#include "Surface.h"

#include <iostream>
#include "learnply_io.h"
#include <Eigen/Dense>



Surface::Surface()
{
}


Surface::~Surface()
{
}

Surface::Surface(FILE *file)
{
	int i, j;
	int elem_count;
	char *elem_name;

	/*** Read in the original PLY object ***/
	static PlyFile *in_ply;
	in_ply = read_ply(file);

	for (i = 0; i < in_ply->num_elem_types; i++) {

		/* prepare to read the i'th list of elements */
		elem_name = setup_element_read_ply(in_ply, i, &elem_count);

		if (equal_strings("vertex", elem_name)) {

			/* create a vertex list to hold all the vertices */
			max_verts = elem_count;

			/* set up for getting vertex elements */

			setup_property_ply(in_ply, &vert_props[0]);
			setup_property_ply(in_ply, &vert_props[1]);
			setup_property_ply(in_ply, &vert_props[2]);
			vert_other = get_other_properties_ply(in_ply,
				offsetof(Vertex_io, other_props));

			/* grab all the vertex elements */
			for (j = 0; j < vertices.size(); j++) {
				Vertex_io vert;
				get_element_ply(in_ply, (void *)&vert);

				/* copy info from the "vert" structure */
				vertices[j] = new Vertex(vert.x, vert.y, vert.z);
				vertices[j]->other_props = vert.other_props;
			}
		}
		else if (equal_strings("face", elem_name)) {

			/* create a list to hold all the face elements */
			max_tris = elem_count;

			/* set up for getting face elements */
			setup_property_ply(in_ply, &face_props[0]);
			face_other = get_other_properties_ply(in_ply, offsetof(Face_io, other_props));

			/* grab all the face elements */
			for (j = 0; j < elem_count; j++) {
				Face_io face;
				get_element_ply(in_ply, (void *)&face);

				if (face.nverts != 3) {
					fprintf(stderr, "Face has %d vertices (should be three).\n",
						face.nverts);
					exit(-1);
				}

				/* copy info from the "face" structure */
				triangles[j] = new Triangle;
				triangles[j]->nverts = 3;
				triangles[j]->verts[0] = (Vertex *)face.verts[0];
				triangles[j]->verts[1] = (Vertex *)face.verts[1];
				triangles[j]->verts[2] = (Vertex *)face.verts[2];
				triangles[j]->other_props = face.other_props;

			}
		}
		else
			get_other_element_ply(in_ply);
	}

	/* close the file */
	close_ply(in_ply);

	/* fix up vertex pointers in triangles */
	for (i = 0; i < triangles.size(); i++) {
		triangles[i]->verts[0] = vertices[(int)(size_t)triangles[i]->verts[0]];
		triangles[i]->verts[1] = vertices[(int)(size_t)triangles[i]->verts[1]];
		triangles[i]->verts[2] = vertices[(int)(size_t)triangles[i]->verts[2]];
	}

	/* get rid of triangles that use the same vertex more than once */

	for (i = triangles.size() - 1; i >= 0; i--) {

		Triangle *tri = triangles[i];
		Vertex *v0 = tri->verts[0];
		Vertex *v1 = tri->verts[1];
		Vertex *v2 = tri->verts[2];

		if (v0 == v1 || v1 == v2 || v2 == v0) {
			triangles.erase(triangles.begin() + i);
			triangles[i] = triangles[triangles.size()];
		}
	}
}

void Surface::initialize() {
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	seed = -1;
}

void Surface::finalize() {
	int i;

	for (i = 0; i<triangles.size(); i++) {
		free(triangles[i]->other_props);
		free(triangles[i]);
	}
	for (i = 0; i<edges.size(); i++) {
		free(edges[i]->tris);
		free(edges[i]);
	}
	for (i = 0; i<vertices.size(); i++) {
		free(vertices[i]->tris);
		free(vertices[i]->other_props);
		free(vertices[i]);
	}

	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
f1    - face that we're looking to share with
v1,v2 - two vertices of f1 that define edge

Exit:
return the matching face, or NULL if there is no such face
******************************************************************************/

Triangle *Surface::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
	int i, j;
	Triangle *f2;
	Triangle *adjacent = NULL;

	/* look through all faces of the first vertex */

	for (i = 0; i < v1->ntris; i++) {
		f2 = v1->tris[i];
		if (f2 == f1)
			continue;
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < f2->nverts; j++) {

			/* look for a match */
			if (f2->verts[j] == v2) {

#if 0
				/* watch out for triple edges */

				if (adjacent != NULL) {

					fprintf(stderr, "model has triple edges\n");

					fprintf(stderr, "face 1: ");
					for (k = 0; k < f1->nverts; k++)
						fprintf(stderr, "%d ", f1->iverts[k]);
					fprintf(stderr, "\nface 2: ");
					for (k = 0; k < f2->nverts; k++)
						fprintf(stderr, "%d ", f2->iverts[k]);
					fprintf(stderr, "\nface 3: ");
					for (k = 0; k < adjacent->nverts; k++)
						fprintf(stderr, "%d ", adjacent->iverts[k]);
					fprintf(stderr, "\n");

				}

				/* if we've got a match, remember this face */
				adjacent = f2;
#endif

#if 1
				/* if we've got a match, return this face */
				return (f2);
#endif

			}
		}
	}

	return (adjacent);
}


/******************************************************************************
Create an edge.

Entry:
v1,v2 - two vertices of f1 that define edge
******************************************************************************/

void Surface::create_edge(Vertex *v1, Vertex *v2)
{
	int i, j;
	Triangle *f;

	/* create the edge */

	edges[edges.size()] = new Edge;
	Edge *e = edges[edges.size()];
	e->index = edges.size();
	e->verts[0] = v1;
	e->verts[1] = v2;
	e->middle = NULL;

	/* count all triangles that will share the edge, and do this */
	/* by looking through all faces of the first vertex */

	e->ntris = 0;

	for (i = 0; i < v1->ntris; i++) {
		f = v1->tris[i];
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 3; j++) {
			/* look for a match */
			if (f->verts[j] == v2) {
				e->ntris++;
				break;
			}
		}
	}

	/* make room for the face pointers (at least two) */
	if (e->ntris < 2)
		e->tris = new Triangle *[2];
	else
		e->tris = new Triangle *[e->ntris];

	/* create pointers from edges to faces and vice-versa */

	e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

	for (i = 0; i < v1->ntris; i++) {

		f = v1->tris[i];

		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 3; j++)
			if (f->verts[j] == v2) {

				e->tris[e->ntris] = f;
				e->ntris++;

				if (f->verts[(j + 1) % 3] == v1)
					f->edges[j] = e;
				else if (f->verts[(j + 2) % 3] == v1)
					f->edges[(j + 2) % 3] = e;
				else {
					fprintf(stderr, "Non-recoverable inconsistancy in create_edge()\n");
					exit(-1);
				}

				break;  /* we'll only find one instance of v2 */
			}

	}
}


/******************************************************************************
Create edges.
******************************************************************************/

void Surface::create_edges()
{
	int i, j;
	Triangle *f;
	Vertex *v1, *v2;
	double count = 0;

	/* count up how many edges we may require */

	for (i = 0; i < triangles.size(); i++) {
		f = triangles[i];
		for (j = 0; j < f->nverts; j++) {
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			Triangle *result = find_common_edge(f, v1, v2);
			if (result)
				count += 0.5;
			else
				count += 1;
		}
	}

	/*
	printf ("counted %f edges\n", count);
	*/

	/* create space for edge list */

	max_edges = (int)(count + 10);  /* leave some room for expansion */

	/* zero out all the pointers from faces to edges */

	for (i = 0; i < triangles.size(); i++)
		for (j = 0; j < 3; j++)
			triangles[i]->edges[j] = NULL;

	/* create all the edges by examining all the triangles */

	for (i = 0; i < triangles.size(); i++) {
		f = triangles[i];
		for (j = 0; j < 3; j++) {
			/* skip over edges that we've already created */
			if (f->edges[j])
				continue;
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			create_edge(v1, v2);
		}
	}
}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Surface::vertex_to_tri_ptrs()
{
	int i, j;
	Triangle *f;
	Vertex *v;

	/* zero the count of number of pointers to faces */

	for (i = 0; i < vertices.size(); i++)
		vertices[i]->max_tris = 0;

	/* first just count all the face pointers needed for each vertex */

	for (i = 0; i < triangles.size(); i++) {
		f = triangles[i];
		for (j = 0; j < f->nverts; j++)
			f->verts[j]->max_tris++;
	}

	/* allocate memory for face pointers of vertices */

	for (i = 0; i < vertices.size(); i++) {
		vertices[i]->tris = (Triangle **)
			malloc(sizeof(Triangle *) * vertices[i]->max_tris);
		vertices[i]->ntris = 0;
	}

	/* now actually create the face pointers */

	for (i = 0; i < triangles.size(); i++) {
		f = triangles[i];
		for (j = 0; j < f->nverts; j++) {
			v = f->verts[j];
			v->tris[v->ntris] = f;
			v->ntris++;
		}
	}
}


/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/

Triangle *Surface::other_triangle(Edge *edge, Triangle *tri)
{
	/* search for any other triangle */

	for (int i = 0; i < edge->ntris; i++)
		if (edge->tris[i] != tri)
			return (edge->tris[i]);

	/* there is no such other triangle if we get here */
	return (NULL);
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
v - vertex whose face list is to be ordered
******************************************************************************/

void Surface::order_vertex_to_tri_ptrs(Vertex *v)
{
	int i, j;
	Triangle *f;
	Triangle *fnext;
	int nf;
	int vindex;
	int boundary;
	int count;

	nf = v->ntris;
	f = v->tris[0];

	/* go backwards (clockwise) around faces that surround a vertex */
	/* to find out if we reach a boundary */

	boundary = 0;

	for (i = 1; i <= nf; i++) {

		/* find reference to v in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[j] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #1\n");
			exit(-1);
		}

		/* corresponding face is the previous one around v */
		fnext = other_triangle(f->edges[vindex], f);

		/* see if we've reached a boundary, and if so then place the */
		/* current face in the first position of the vertice's face list */

		if (fnext == NULL) {
			/* find reference to f in v */
			for (j = 0; j < v->ntris; j++)
				if (v->tris[j] == f) {
					v->tris[j] = v->tris[0];
					v->tris[0] = f;
					break;
				}
			boundary = 1;
			break;
		}

		f = fnext;
	}

	/* now walk around the faces in the forward direction and place */
	/* them in order */

	f = v->tris[0];
	count = 0;

	for (i = 1; i < nf; i++) {

		/* find reference to vertex in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[(j + 1) % f->nverts] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #2\n");
			exit(-1);
		}

		/* corresponding face is next one around v */
		fnext = other_triangle(f->edges[vindex], f);

		/* break out of loop if we've reached a boundary */
		count = i;
		if (fnext == NULL) {
			break;
		}

		/* swap the next face into its proper place in the face list */
		for (j = 0; j < v->ntris; j++)
			if (v->tris[j] == fnext) {
				v->tris[j] = v->tris[i];
				v->tris[i] = fnext;
				break;
			}

		f = fnext;
	}
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
f - face whose vertex list is to be searched
v - vertex to return reference to

Exit:
returns index in face's list, or -1 if vertex not found
******************************************************************************/

int Surface::face_to_vertex_ref(Triangle *f, Vertex *v)
{
	int j;
	int vindex = -1;

	for (j = 0; j < f->nverts; j++)
		if (f->verts[j] == v) {
			vindex = j;
			break;
		}

	return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/

void Surface::create_pointers()
{
	int i;

	/* index the vertices and triangles */

	for (i = 0; i < vertices.size(); i++)
		vertices[i]->index = i;

	for (i = 0; i < triangles.size(); i++)
		triangles[i]->index = i;

	/* create pointers from vertices to triangles */
	vertex_to_tri_ptrs();

	/* make edges */
	create_edges();


	/* order the pointers from vertices to faces */
	for (i = 0; i < vertices.size(); i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "ordering %d of %d vertices\n", i, vertices.size());
		order_vertex_to_tri_ptrs(vertices[i]);

	}
	/* index the edges */

	for (i = 0; i < edges.size(); i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
		edges[i]->index = i;
	}

}

void Surface::calc_bounding_sphere()
{
	unsigned int i;

	for (i = 0; i<vertices.size(); i++) {
		if (i == 0) {
			min.set(vertices[i]->x, vertices[i]->y, vertices[i]->z);
			max.set(vertices[i]->x, vertices[i]->y, vertices[i]->z);
		}
		else {
			if (vertices[i]->x < min.entry[0])
				min.entry[0] = vertices[i]->x;
			if (vertices[i]->x > max.entry[0])
				max.entry[0] = vertices[i]->x;
			if (vertices[i]->y < min.entry[1])
				min.entry[1] = vertices[i]->y;
			if (vertices[i]->y > max.entry[1])
				max.entry[1] = vertices[i]->y;
			if (vertices[i]->z < min.entry[2])
				min.entry[2] = vertices[i]->z;
			if (vertices[i]->z > max.entry[2])
				max.entry[2] = vertices[i]->z;
		}
	}
	center = (min + max) * 0.5;
	radius = length(center - min);
}

void Surface::calc_moments() {
	double sum_vx = 0.0, sum_vy = 0.0, sum_vz = 0.0;
	for (int i = 0; i<vertices.size(); i++) {
		sum_vx += vertices[i]->x;
		sum_vy += vertices[i]->y;
		sum_vz += vertices[i]->z;
	}
	sum_vx /= vertices.size();
	sum_vy /= vertices.size();
	sum_vz /= vertices.size();
	gravity_center = icVector3(sum_vx, sum_vy, sum_vz);

	Eigen::Matrix3d variance = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();

	for (int i = 0; i<vertices.size(); i++) {
		double x_v = vertices[i]->x - gravity_center.x, y_v = vertices[i]->y - gravity_center.y, z_v = vertices[i]->z - gravity_center.z;
		double x_n = vertices[i]->normal.x - gravity_center.x, y_n = vertices[i]->normal.y - gravity_center.y, z_n = vertices[i]->normal.z - gravity_center.z;
		Eigen::Matrix3d temp_v;
		temp_v << x_v * x_v, x_v * y_v, x_v * z_v, y_v * x_v, y_v * y_v, y_v * z_v, z_v * x_v, z_v * y_v, z_v * z_v;
		Eigen::Matrix3d temp_n;
		temp_n << x_n * x_n, x_n * y_n, x_n * z_n, y_n * x_n, y_n * y_n, y_n * z_n, z_n * x_n, z_n * y_n, z_n * z_n;
		variance += temp_v;
		normal += temp_n;
	}
	Eigen::EigenSolver<Eigen::MatrixXd> es_v(variance);
	Eigen::EigenSolver<Eigen::MatrixXd> es_n(normal);

	std::vector <double> eigen_v(es_v.eigenvectors().size());
	Eigen::Map<Eigen::MatrixXd>(eigen_v.data(), es_v.eigenvectors().rows(), es_v.eigenvectors().cols()) = es_v.eigenvectors().real();

	std::vector <double> eigen_n(es_n.eigenvectors().size());
	Eigen::Map<Eigen::MatrixXd>(eigen_n.data(), es_n.eigenvectors().rows(), es_n.eigenvectors().cols()) = es_n.eigenvectors().real();

	double d_x_p = 0.0;
	double d_x_n = 0.0;
	double d_y_p = 0.0;
	double d_y_n = 0.0;
	double d_z_p = 0.0;
	double d_z_n = 0.0;

	for (int i = 0; i<vertices.size(); i++) {
		icVector3 vertex = icVector3(vertices[i]->x, vertices[i]->y, vertices[i]->z);
		double ev_x = vertex.x * eigen_v[0] + vertex.y * eigen_v[1] + vertex.z * eigen_v[2];
		double ev_y = vertex.x * eigen_v[3] + vertex.y * eigen_v[4] + vertex.z * eigen_v[5];
		double ev_z = vertex.x * eigen_v[6] + vertex.y * eigen_v[7] + vertex.z * eigen_v[8];
		if (ev_x > 0.0) {
			d_x_p = std::max(d_x_p, ev_x);
		}
		else {
			d_x_n = std::min(d_x_n, ev_x);
		}
		if (ev_y > 0.0) {
			d_y_p = std::max(d_y_p, ev_y);
		}
		else {
			d_y_n = std::min(d_y_n, ev_y);
		}
		if (ev_z > 0.0) {
			d_z_p = std::max(d_z_p, ev_z);
		}
		else {
			d_z_n = std::min(d_z_n, ev_z);
		}
	}

	moment_vertices = { icVector3(gravity_center.x + d_x_p, gravity_center.y + d_y_p, gravity_center.z + d_z_p),
		icVector3(gravity_center.x + d_x_p, gravity_center.y + d_y_p, gravity_center.z + d_z_n),
		icVector3(gravity_center.x + d_x_p, gravity_center.y + d_y_n, gravity_center.z + d_z_p),
		icVector3(gravity_center.x + d_x_p, gravity_center.y + d_y_n, gravity_center.z + d_z_n),
		icVector3(gravity_center.x + d_x_n, gravity_center.y + d_y_p, gravity_center.z + d_z_p),
		icVector3(gravity_center.x + d_x_n, gravity_center.y + d_y_p, gravity_center.z + d_z_n),
		icVector3(gravity_center.x + d_x_n, gravity_center.y + d_y_n, gravity_center.z + d_z_p),
		icVector3(gravity_center.x + d_x_n, gravity_center.y + d_y_n, gravity_center.z + d_z_n)
	};

	double n_x_p = 0.0;
	double n_x_n = 0.0;
	double n_y_p = 0.0;
	double n_y_n = 0.0;
	double n_z_p = 0.0;
	double n_z_n = 0.0;

	for (int i = 0; i<vertices.size(); i++) {
		icVector3 vertex = icVector3(vertices[i]->x, vertices[i]->y, vertices[i]->z);
		double ev_x = vertex.x * eigen_n[0] + vertex.y * eigen_n[1] + vertex.z * eigen_n[2];
		double ev_y = vertex.x * eigen_n[3] + vertex.y * eigen_n[4] + vertex.z * eigen_n[5];
		double ev_z = vertex.x * eigen_n[6] + vertex.y * eigen_n[7] + vertex.z * eigen_n[8];
		if (ev_x > 0.0) {
			n_x_p = std::max(n_x_p, ev_x);
		}
		else {
			n_x_n = std::min(n_x_n, ev_x);
		}
		if (ev_y > 0.0) {
			n_y_p = std::max(n_y_p, ev_y);
		}
		else {
			n_y_n = std::min(n_y_n, ev_y);
		}
		if (ev_z > 0.0) {
			n_z_p = std::max(n_z_p, ev_z);
		}
		else {
			n_z_n = std::min(n_z_n, ev_z);
		}
	}

	normal_vertices = { icVector3(gravity_center.x + n_x_p, gravity_center.y + n_y_p, gravity_center.z + n_z_p),
		icVector3(gravity_center.x + n_x_p, gravity_center.y + n_y_p, gravity_center.z + n_z_n),
		icVector3(gravity_center.x + n_x_p, gravity_center.y + n_y_n, gravity_center.z + n_z_p),
		icVector3(gravity_center.x + n_x_p, gravity_center.y + n_y_n, gravity_center.z + n_z_n),
		icVector3(gravity_center.x + n_x_n, gravity_center.y + n_y_p, gravity_center.z + n_z_p),
		icVector3(gravity_center.x + n_x_n, gravity_center.y + n_y_p, gravity_center.z + n_z_n),
		icVector3(gravity_center.x + n_x_n, gravity_center.y + n_y_n, gravity_center.z + n_z_p),
		icVector3(gravity_center.x + n_x_n, gravity_center.y + n_y_n, gravity_center.z + n_z_n)
	};

}

void Surface::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i = 0; i<egdes.size(); i++) {
		v1.set(edges[i]->verts[0]->x, edges[i]->verts[0]->y, edges[i]->verts[0]->z);
		v2.set(edges[i]->verts[1]->x, edges[i]->verts[1]->y, edges[i]->verts[1]->z);
		edges[i]->length = length(v1 - v2);
	}
}

void Surface::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2;
	Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (i = 0; i<triangles.size(); i++) {
		for (j = 0; j<3; j++)
			length[j] = triangles[i]->edges[j]->length;
		double temp_s = (length[0] + length[1] + length[2]) / 2.0;
		triangles[i]->area = sqrt(temp_s*(temp_s - length[0])*(temp_s - length[1])*(temp_s - length[2]));

		area += triangles[i]->area;
		temp_t = triangles[i];
		v1.set(vertices[triangles[i]->verts[0]->index]->x, vertices[triangles[i]->verts[0]->index]->y, vertices[triangles[i]->verts[0]->index]->z);
		v2.set(vertices[triangles[i]->verts[1]->index]->x, vertices[triangles[i]->verts[1]->index]->y, vertices[triangles[i]->verts[1]->index]->z);
		v0.set(vertices[triangles[i]->verts[2]->index]->x, vertices[triangles[i]->verts[2]->index]->y, vertices[triangles[i]->verts[2]->index]->z);
		triangles[i]->normal = cross(v0 - v1, v2 - v1);
		normalize(triangles[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i = 0; i<triangles.size(); i++) {
		icVector3 cent(vertices[triangles[i]->verts[0]->index]->x, vertices[triangles[i]->verts[0]->index]->y, vertices[triangles[i]->verts[0]->index]->z);
		signedvolume += dot(test - cent, triangles[i]->normal)*triangles[i]->area;
	}
	signedvolume /= area;
	if (signedvolume<0)
		orientation = 0;
	else {
		orientation = 1;
		for (i = 0; i<ntris; i++)
			triangles[i]->normal *= -1.0;
	}
}

bool sortcol(const std::vector<int>& v1, const std::vector<int>& v2) {
	return ((v1[1] < v2[1]) || ((v1[1] == v2[1]) && (v1[2]<v2[2])));
}

void Surface::calc_corner_table()
{
	/* Build corner table */

	int ncorners = 0;

	for (int i = 0; i < triangles.size; i++) {
		for (int k = 0; k < 3; k++) {
			Vertex* v0 = triangles[i]->verts[k];
			Vertex* v1 = triangles[i]->verts[(k + 1) % 3];
			Vertex* v2 = triangles[i]->verts[(k + 2) % 3];
			corners[ncorners] = new Corner;
			corners[ncorners]->index = ncorners;
			corners[ncorners]->triangle = triangles[i];
			corners[ncorners]->vertex = triangles[i]->verts[k];
			corners[ncorners]->edge = triangles[i]->edges[(k + 2) % 3];
			double vec1[] = { v1->x - v0->x, v1->y - v0->y, v1->z - v0->z };
			double vec2[] = { v2->x - v0->x, v2->y - v0->y, v2->z - v0->z };
			double vec1mag = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
			double vec1_norm[] = { vec1[0] / vec1mag, vec1[1] / vec1mag, vec1[2] / vec1mag };
			double vec2mag = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);
			double vec2_norm[] = { vec2[0] / vec2mag, vec2[1] / vec2mag, vec2[2] / vec2mag };
			double res = vec1_norm[0] * vec2_norm[0] + vec1_norm[1] * vec2_norm[1] + vec1_norm[2] * vec2_norm[2];
			corners[ncorners]->angle = acos(res);
			corners[ncorners]->o = NULL;
			if (k > 0) {
				corners[ncorners]->p = corners[ncorners - 1];
				corners[ncorners - 1]->n = corners[ncorners];
			}
			ncorners++;
		}
		corners[ncorners - 3]->p = corners[ncorners - 1];
		corners[ncorners - 1]->n = corners[ncorners - 3];
	}


	std::vector<std::vector<int>> opposite_table;

	for (int i = 0; i < ncorners; i++) {
		Corner* c = corners[i];
		std::vector<int> min_max = { c->index, std::min(c->p->vertex->index, c->n->vertex->index), std::max(c->p->vertex->index, c->n->vertex->index) };
		opposite_table.push_back(min_max);
	}
	sort(opposite_table.begin(), opposite_table.end(), sortcol);
	for (int i = 0; i < opposite_table.size() - 1; i++) {
		if (opposite_table[i][1] == opposite_table[i + 1][1] && opposite_table[i][2] == opposite_table[i + 1][2]) {
			corners[opposite_table[i][0]]->o = corners[opposite_table[i + 1][0]];
			corners[opposite_table[i + 1][0]]->o = corners[opposite_table[i][0]];
		}
	}
}

void Surface::test_corner_table()
{
	bool error = false;
	for (int i = 0; i < corners.size(); i++)
	{
		if (corners[i]->n->p != corners[i]) {
			std::cout << "Next corners previous is not current corner" << std::endl;
			error = true;
		}
		if (corners[i]->n->n->n != corners[i]) {
			std::cout << "Next corner x3 is not current corner" << std::endl;
			error = true;
		}
		if (corners[i]->o != NULL && corners[i]->o->o != corners[i]) {
			std::cout << "Opposite of opposite is not the same" << std::endl;
			error = true;
		}
	}
	if (!error)
		std::cout << "Corner table successfully created. Number of corners: " + std::to_string(corners.size()) << std::endl;
}
