// CS552 Final Project.cpp

#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <gl/GLUT.h>

#include "Surface.h"

#define SCREEN_WIDTH	800
#define SCREEN_HEIGHT	600
#define WINDOW_TITLE	"CS552 Final Project"

Surface* surface;

void display() {
	glEnable(GL_LIGHTING);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();


	glFlush();
	glutSwapBuffers();
}

void timer(int value) {
	glutPostRedisplay();
	glutTimerFunc(1000, timer, 1);
}

void idle()
{
	glutPostRedisplay();
}

void reshape(int w, int h) {

}

void onKeyPress(unsigned char key, int x, int y) {

}

void onMouseButton(int button, int state, int x, int y) {

}

void onMouseDrag(int x, int y) {

}

int main(int argc, char *argv[])
{
	FILE* this_file = fopen("tempmodels/bunny.ply", "r");
	surface = new Surface(this_file);
	fclose(this_file);

	surface->initialize(); // initialize everything
	surface->calc_bounding_sphere();
	surface->calc_face_normals_and_area();
	surface->calc_moments();

	surface->calc_corner_table();
	surface->test_corner_table();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(WINDOW_TITLE);

	glutDisplayFunc(display);
	glutTimerFunc(1000, timer, 0);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(onKeyPress);
	glutMouseFunc(onMouseButton);
	glutMotionFunc(onMouseDrag);
	glutIdleFunc(idle);

	if (glewInit() != GLEW_OK)
	{
		fprintf(stdout, "Failed to initialize GLEW");
		return EXIT_FAILURE;
	}
	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);

	GLfloat global_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glLightModelfv(GL_AMBIENT_AND_DIFFUSE, global_ambient);

	GLfloat aspect = (GLfloat)SCREEN_WIDTH / (GLfloat)SCREEN_HEIGHT;
	glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, aspect, 0.1f, 1000.0f);

	glutMainLoop();

	return EXIT_SUCCESS;
}

