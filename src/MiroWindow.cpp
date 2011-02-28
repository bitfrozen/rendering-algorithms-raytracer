#include "MiroWindow.h"
#include "OpenGL.h"
#include "Miro.h"
#include "Camera.h"
#include "Image.h"
#include <stdlib.h>
#include <time.h>
#include "Scene.h"

#define ANGFACT     1.0
#define LEFT        4
#define MIDDLE      2
#define RIGHT       1

#ifdef WIN32
// disable useless warnings
#pragma warning(disable:4996)
#endif

namespace
{

// Non-member functions used as proxy callbacks to our real C++ member functions
MiroWindow *g_miroWindow;
void display() {g_miroWindow->display();}
void resize(int x,int y) {g_miroWindow->reshape(x,y);}
void keyboard(unsigned char key, int x, int y) {g_miroWindow->keyboard(key,x,y);}
void specialKeys(int key, int x, int y) {g_miroWindow->specialKeys(key,x,y);}
void mouse(int btn,int state,int x,int y) {g_miroWindow->mouse(btn,state,x,y);}
void motion(int x, int y) {g_miroWindow->motion(x,y);}

} // namespace

MiroWindow::MiroWindow(int * argc, char* argv[]) :
    m_scaleFact(0.1f),
	m_microScaleFact(0.01f),
    m_activeButton(0),
    m_mouseX(0),
    m_mouseY(0)
{
	theParam = SCALE_FACT;
    // Initialize GLUT
    glutInit(argc, argv);

    // Create the window
    glutInitWindowSize(g_image->width(), g_image->height());
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowPosition(200, 200);
    glutCreateWindow("miro");

    // Initialize some OpenGL state
    glClearColor(0.25f, 0.25f, 0.25f, 1);
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glShadeModel(GL_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // draw outlines only
}

void
MiroWindow::mainLoop()
{
    // Setup callback functions
    g_miroWindow = this;
    glutDisplayFunc(::display);
    glutKeyboardFunc(::keyboard);
    glutMouseFunc(::mouse);
    glutMotionFunc(::motion);
    glutReshapeFunc(::resize);
	glutSpecialFunc(::specialKeys);
    
    // Start the glut main loop, never returns
    glutMainLoop();
}


void
MiroWindow::display()
{
    g_camera->click(g_scene, g_image); // take a snapshot of the scene

    glFinish(); // flush the openGL pipeline
}


void
MiroWindow::motion(int x, int y)
{
    int dx, dy;     // change in mouse coordinates

    dx = x - m_mouseX;     // change in mouse coords
    dy = y - m_mouseY;

    if (m_activeButton & LEFT)
    {
        float xfact = -ANGFACT*dy;
        float yfact = -ANGFACT*dx;
        // construct a coordinate system from up and viewdir
        Vector3 vRight = cross(g_camera->viewDir(), g_camera->up());
        // now rotate everything
        Vector3 v = g_camera->viewDir();
        v.rotate(xfact*PI/180., vRight);
        v.rotate(yfact*PI/180., g_camera->up());
        g_camera->setViewDir(v);
    }

    m_mouseX = x;          // new current position
    m_mouseY = y;

    glutPostRedisplay();
}


void
MiroWindow::mouse(int button, int state, int x, int y)
{
    int b; // LEFT, MIDDLE, or RIGHT

    switch (button)
    {
        case GLUT_LEFT_BUTTON:
            b = LEFT;
        break;
        
        case GLUT_MIDDLE_BUTTON:
            b = MIDDLE;
        break;
        
        case GLUT_RIGHT_BUTTON:
            b = RIGHT;
        break;
        
        default:
            b = 0;
    }

    if (state == GLUT_DOWN)
    {
        m_mouseX = x;
        m_mouseY = y;
        m_activeButton |= b;       /* set the proper bit   */
    }
    else
        m_activeButton &= ~b;      /* clear the proper bit */
}

void MiroWindow::increaseParam()
{
	switch (theParam)
	{
	case SCALE_FACT:
		{
			m_scaleFact *= 2.0f;
			m_microScaleFact = m_scaleFact * 0.1f;
			printf("Scale Factor (+/-/up/dn): %.6f\r", m_scaleFact);
			fflush(stdout);
			break;
		}
	case FOV:
		{
			g_camera->setFOV(std::min(179.f, g_camera->FOV() + m_scaleFact));
			printf("Camera - (F)OV: %.3f\r", g_camera->FOV());
			fflush(stdout);
			break;
		}
	case FOCUS_PLANE:
		{
			g_camera->setFocusPlane(g_camera->focusPlane() + m_scaleFact);
			printf("Camera - F(o)cus Distance: %.3f\r", g_camera->focusPlane());
			fflush(stdout);
			break;
		}
	case APERTURE:
		{
			g_camera->setAperture(g_camera->aperture() + m_scaleFact);
			printf("Camera - A(p)erture: %.3f\r", g_camera->aperture());
			fflush(stdout);
			break;
		}
	case NUM_PATHS:
		{
			g_scene->setNumPaths(g_scene->numPaths() + 1);
			printf("Scene  - Number of Pat(h)s: %i\r", g_scene->numPaths());
			fflush(stdout);
			break;
		}
	case NUM_BOUNCES:
		{
			g_scene->setMaxBounces(g_scene->maxBounces() + 1);
			printf("Scene  - Number of (B)ounces: %i\r", g_scene->maxBounces());
			fflush(stdout);
			break;
		}
	case MIN_SUBDIVS:
		{
			g_scene->setMinSubdivs(std::min(g_scene->maxSubdivs(), g_scene->minSubdivs() + 1));
			printf("Scene  - Minimum S(u)bdivs: %i\r", g_scene->minSubdivs());
			fflush(stdout);
			break;
		}
	case MAX_SUBDIVS:
		{
			g_scene->setMaxSubdivs(g_scene->maxSubdivs() + 1);
			printf("Scene  - Maximum Subdi(v)s: %i\r", g_scene->maxSubdivs());
			fflush(stdout);
			break;
		}
	case NOISE_THRESHOLD:
		{
			g_scene->setNoise(std::min(0.999f, g_scene->noise() + m_scaleFact));
			printf("Scene  - (N)oise Threshold: %.6f\r", g_scene->noise());
			fflush(stdout);
			break;
		}
	case SHUTTER_SPEED:
		{
			g_camera->setShutterSpeed(std::min(1.0f, g_camera->shutterSpeed() + m_scaleFact));
			printf("Camera - Shutter sp(e)ed: %.3ff\r", g_camera->shutterSpeed());
			fflush(stdout);
			break;
		}
	}
}

void MiroWindow::decreaseParam()
{
	switch (theParam)
	{
	case SCALE_FACT:
		{
			m_scaleFact /= 2.0f;
			m_microScaleFact = m_scaleFact * 0.1f;
			printf("Scale Factor (+/-/up/dn): %.6f\r", m_scaleFact);
			fflush(stdout);
			break;
		}
	case FOV:
		{
			g_camera->setFOV(std::max(0.5f, g_camera->FOV() - m_scaleFact));
			printf("Camera - (F)OV: %.3f\r", g_camera->FOV());
			fflush(stdout);
			break;
		}
	case FOCUS_PLANE:
		{
			g_camera->setFocusPlane(std::max(0.001f, g_camera->focusPlane() - m_scaleFact));
			printf("Camera - F(o)cus Distance: %.3f\r", g_camera->focusPlane());
			fflush(stdout);
			break;
		}
	case APERTURE:
		{
			g_camera->setAperture(std::max(0.0f, g_camera->aperture() - m_scaleFact));
			printf("Camera - A(p)erture: %.3f\r", g_camera->aperture());
			fflush(stdout);
			break;
		}
	case NUM_PATHS:
		{
			g_scene->setNumPaths(std::max(1, g_scene->numPaths() - 1));
			printf("Scene  - Number of Pat(h)s: %i\r", g_scene->numPaths());
			fflush(stdout);
			break;
		}
	case NUM_BOUNCES:
		{
			g_scene->setMaxBounces(std::max(1, g_scene->maxBounces() - 1));
			printf("Scene  - Number of (B)ounces: %i\r", g_scene->maxBounces());
			fflush(stdout);
			break;
		}
	case MIN_SUBDIVS:
		{
			g_scene->setMinSubdivs(std::max(1, g_scene->minSubdivs() - 1));
			printf("Scene  - Minimum S(u)bdivs: %i\r", g_scene->minSubdivs());
			fflush(stdout);
			break;
		}
	case MAX_SUBDIVS:
		{
			g_scene->setMaxSubdivs(std::max(g_scene->minSubdivs(), g_scene->maxSubdivs() - 1));
			printf("Scene  - Maximum Subdi(v)s: %i\r", g_scene->maxSubdivs());
			fflush(stdout);
			break;
		}
	case NOISE_THRESHOLD:
		{
			g_scene->setNoise(std::max(0.0f, g_scene->noise() - m_scaleFact));
			printf("Scene  - (N)oise Threshold: %.6f\r", g_scene->noise());
			fflush(stdout);
			break;
		}
	case SHUTTER_SPEED:
		{
			g_camera->setShutterSpeed(std::max(0.0f, g_camera->shutterSpeed() - m_scaleFact));
			printf("Camera - Shutter sp(e)ed: %.3ff\r", g_camera->shutterSpeed());
			fflush(stdout);
			break;
		}
	}
}

void MiroWindow::specialKeys(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:
		{
			switch (theParam)
			{
			case SCALE_FACT:
				{
					m_scaleFact += m_microScaleFact;
					printf("Scale Factor (+/-/up/dn): %.6f\r", m_scaleFact);
					fflush(stdout);
					break;
				}
			case FOV:
				{
					g_camera->setFOV(std::max(0.5f, g_camera->FOV() + m_microScaleFact));
					printf("Camera - (F)OV: %.3f\r", g_camera->FOV());
					fflush(stdout);
					break;
				}
			case FOCUS_PLANE:
				{
					g_camera->setFocusPlane(std::max(0.001f, g_camera->focusPlane() + m_microScaleFact));
					printf("Camera - F(o)cus Distance: %.3f\r", g_camera->focusPlane());
					fflush(stdout);
					break;
				}
			case APERTURE:
				{
					g_camera->setAperture(std::max(0.0f, g_camera->aperture() + m_microScaleFact));
					printf("Camera - A(p)erture: %.3f\r", g_camera->aperture());
					fflush(stdout);
					break;
				}
			case NUM_PATHS:
				{
					g_scene->setNumPaths(std::max(1, g_scene->numPaths() + 1));
					printf("Scene  - Number of Pat(h)s: %i\r", g_scene->numPaths());
					fflush(stdout);
					break;
				}
			case NUM_BOUNCES:
				{
					g_scene->setMaxBounces(std::max(1, g_scene->maxBounces() + 1));
					printf("Scene  - Number of (B)ounces: %i\r", g_scene->maxBounces());
					fflush(stdout);
					break;
				}
			case MIN_SUBDIVS:
				{
					g_scene->setMinSubdivs(std::max(1, g_scene->minSubdivs() + 1));
					printf("Scene  - Minimum S(u)bdivs: %i\r", g_scene->minSubdivs());
					fflush(stdout);
					break;
				}
			case MAX_SUBDIVS:
				{
					g_scene->setMaxSubdivs(std::max(g_scene->minSubdivs(), g_scene->maxSubdivs() + 1));
					printf("Scene  - Maximum Subdi(v)s: %i\r", g_scene->maxSubdivs());
					fflush(stdout);
					break;
				}
			case NOISE_THRESHOLD:
				{
					g_scene->setNoise(std::max(0.0f, g_scene->noise() + m_microScaleFact));
					printf("Scene  - (N)oise Threshold: %.6f\r", g_scene->noise());
					fflush(stdout);
					break;
				}
			case SHUTTER_SPEED:
				{
					g_camera->setShutterSpeed(std::min(1.0f, g_camera->shutterSpeed() + m_microScaleFact));
					printf("Camera - Shutter sp(e)ed: %.3ff\r", g_camera->shutterSpeed());
					fflush(stdout);
					break;
				}
			}
			break;
		}
	case GLUT_KEY_DOWN:
		{
			switch (theParam)
			{
			case SCALE_FACT:
				{
					m_scaleFact -= m_microScaleFact;
					m_scaleFact = max(0.001f, m_scaleFact);
					printf("Scale Factor (+/-/up/dn): %.6f\r", m_scaleFact);
					fflush(stdout);
					break;
				}
			case FOV:
				{
					g_camera->setFOV(std::max(0.5f, g_camera->FOV() - m_microScaleFact));
					printf("Camera - (F)OV: %.3f\r", g_camera->FOV());
					fflush(stdout);
					break;
				}
			case FOCUS_PLANE:
				{
					g_camera->setFocusPlane(std::max(0.001f, g_camera->focusPlane() - m_microScaleFact));
					printf("Camera - F(o)cus Distance: %.3f\r", g_camera->focusPlane());
					fflush(stdout);
					break;
				}
			case APERTURE:
				{
					g_camera->setAperture(std::max(0.0f, g_camera->aperture() - m_microScaleFact));
					printf("Camera - A(p)erture: %.3f\r", g_camera->aperture());
					fflush(stdout);
					break;
				}
			case NUM_PATHS:
				{
					g_scene->setNumPaths(std::max(1, g_scene->numPaths() - 1));
					printf("Scene  - Number of Pat(h)s: %i\r", g_scene->numPaths());
					fflush(stdout);
					break;
				}
			case NUM_BOUNCES:
				{
					g_scene->setMaxBounces(std::max(1, g_scene->maxBounces() - 1));
					printf("Scene  - Number of (B)ounces: %i\r", g_scene->maxBounces());
					fflush(stdout);
					break;
				}
			case MIN_SUBDIVS:
				{
					g_scene->setMinSubdivs(std::max(1, g_scene->minSubdivs() - 1));
					printf("Scene  - Minimum S(u)bdivs: %i\r", g_scene->minSubdivs());
					fflush(stdout);
					break;
				}
			case MAX_SUBDIVS:
				{
					g_scene->setMaxSubdivs(std::max(g_scene->minSubdivs(), g_scene->maxSubdivs() - 1));
					printf("Scene  - Maximum Subdi(v)s: %i\r", g_scene->maxSubdivs());
					fflush(stdout);
					break;
				}
			case NOISE_THRESHOLD:
				{
					g_scene->setNoise(std::max(0.0f, g_scene->noise() - m_microScaleFact));
					printf("Scene  - (N)oise Threshold: %.6f\r", g_scene->noise());
					fflush(stdout);
					break;
				}
			case SHUTTER_SPEED:
				{
					g_camera->setShutterSpeed(std::min(1.0f, g_camera->shutterSpeed() - m_microScaleFact));
					printf("Camera - Shutter sp(e)ed: %.3ff\r", g_camera->shutterSpeed());
					fflush(stdout);
					break;
				}
			}
			break;
		}
	default:
		break;
	}
	glutPostRedisplay();
}

void MiroWindow::keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 'i':
        case 'I':
        {
            char str[1024];
            sprintf(str, "miro_%d.ppm", time(0));
            if (g_camera->isOpenGL())
            {
                unsigned char* buf = new unsigned char[g_image->width()*g_image->height()*3];
                glReadPixels(0, 0, g_image->width(), g_image->height(),
                             GL_RGB, GL_UNSIGNED_BYTE, buf);
                g_image->writePPM(str, buf, g_image->width(), g_image->height());
            }
            else
            {
                g_image->writePPM(str);
            }
            break;
        }

        case 'r':
        case 'R':
            g_camera->setRenderer(Camera::RENDER_RAYTRACE);
        break;

        case 'g':
        case 'G':
            g_camera->setRenderer(Camera::RENDER_OPENGL);
        break;

        case '+':
			increaseParam();            
        break;

        case '-':
            decreaseParam();
        break;

        case 'w':
        case 'W':
            g_camera->setEye(g_camera->eye() + m_scaleFact*g_camera->viewDir());
        break;

        case 's':
        case 'S':
            g_camera->setEye(g_camera->eye() - m_scaleFact*g_camera->viewDir());
        break;

        case 'q':
        case 'Q':
            g_camera->setEye(g_camera->eye() + m_scaleFact*g_camera->up());
        break;

        case 'z':
        case 'Z':
            g_camera->setEye(g_camera->eye() - m_scaleFact*g_camera->up());
        break;

        case 'a':
        case 'A':
        {
            Vector3 vRight = cross(g_camera->viewDir(), g_camera->up());
            g_camera->setEye(g_camera->eye() - m_scaleFact*vRight);
            break;
        }

        case 'd':
        case 'D':
        {
            Vector3 vRight = cross(g_camera->viewDir(), g_camera->up());
            g_camera->setEye(g_camera->eye() + m_scaleFact*vRight);
            break;
        }

		case 'm':
		case 'M':
			{
				printf("Scene settings and parameters:\n");
				printf("------------------------------------------------------------------------\n");
				printf("Camera - Position:\tx: %.3f\ty: %.3f\tz: %.3f\n", g_camera->m_eye.x, g_camera->m_eye.y, g_camera->m_eye.z);
				printf("Camera - ViewDir: \tx: %.3f\ty: %.3f\tz: %.3f\n", g_camera->m_viewDir.x, g_camera->m_viewDir.y, g_camera->m_viewDir.z);
				printf("Camera - (F)OV: %.1f\n", g_camera->FOV());
				printf("Camera - F(o)cus Distance: %.3f\n", g_camera->focusPlane());
				printf("Camera - A(p)erture: %.3f\n", g_camera->aperture());
				printf("Camera - Shutter sp(e)ed: %.3ff\n", g_camera->shutterSpeed());
				char *token;
				if (g_scene->m_pathTrace) token = "Yes";
				else token = "No";			
				printf("Scene  - Pa(t)h Tracing: %s.\n", token);
				printf("Scene  - Number of Pat(h)s: %i\n", g_scene->numPaths());
				printf("Scene  - Number of (B)ounces: %i\n", g_scene->maxBounces());
				printf("Scene  - Minimum S(u)bdivs: %i\n", g_scene->minSubdivs());
				printf("Scene  - Maximum Subdi(v)s: %i\n", g_scene->maxSubdivs());
				printf("Scene  - (N)oise Threshold: %.6f\n", g_scene->noise());
				break;
			}

		case 'e':
		case 'E':
			{
				printf("\n");
				if (theParam == SHUTTER_SPEED)
				{
					theParam = SCALE_FACT;
				}				
				else
				{
					theParam = SHUTTER_SPEED;
					printf("Camera - Shutter sp(e)ed: %.3ff\r", g_camera->shutterSpeed());
					fflush(stdout);
				}
				break;
			}

		case 'f':
		case 'F':
			{
				printf("\n");
				if (theParam == FOV)
				{
					theParam = SCALE_FACT;
				}				
				else
				{
					theParam = FOV;
					printf("Camera - (F)OV: %.1f\r", g_camera->m_fov);
					fflush(stdout);
				}
				break;
			}

		case 'o':
		case 'O':
			{
				printf("\n");
				if (theParam == FOCUS_PLANE)
				{
					theParam = SCALE_FACT;
				} 
				else
				{
					theParam = FOCUS_PLANE;
					printf("Camera - F(o)cus Distance: %.3f\r", g_camera->focusPlane());
					fflush(stdout);
				}
				break;
			}

		case 'p':
		case 'P':
			{
				printf("\n");
				if (theParam == APERTURE)
				{
					theParam = SCALE_FACT;
				} 
				else
				{
					theParam = APERTURE;
					printf("Camera - A(p)erture: %.3f\r", g_camera->aperture());
					fflush(stdout);
				}
				break;
			}

		case 't':
		case 'T':
			{
				printf("\n");
				if (g_scene->m_pathTrace)
				{
					g_scene->m_pathTrace = false;
					char *token;
					if (g_scene->m_pathTrace) token = "Yes";
					else token = "No";			
					printf("Scene  - Pa(t)h Tracing: %s.\r", token);
					fflush(stdout);
				} 
				else
				{
					g_scene->m_pathTrace = true;
					char *token;
					if (g_scene->m_pathTrace) token = "Yes";
					else token = "No";			
					printf("Scene  - Pa(t)h Tracing: %s.\r", token);
					fflush(stdout);
				}
				break;
			}

		case 'h':
		case 'H':
			{
				printf("\n");
				if (theParam == NUM_PATHS)
				{
					theParam = SCALE_FACT;
				}				
				else
				{
					theParam = NUM_PATHS;
					printf("Scene  - Number of Pat(h)s: %i\r", g_scene->numPaths());
					fflush(stdout);
				}
				break;
			}

		case 'b':
		case 'B':
			{
				printf("\n");
				if (theParam == NUM_BOUNCES)
				{
					theParam = SCALE_FACT;
				}				
				else
				{
					theParam = NUM_BOUNCES;
					printf("Scene  - Number of (B)ounces: %i\r", g_scene->maxBounces());
					fflush(stdout);
				}
				break;
			}

		case 'u':
		case 'U':
			{
				printf("\n");
				if (theParam == MIN_SUBDIVS)
				{
					theParam = SCALE_FACT;
				}				
				else
				{
					theParam = MIN_SUBDIVS;
					printf("Scene  - Minimum S(u)bdivs: %i\r", g_scene->minSubdivs());
					fflush(stdout);
				}
				break;
			}

		case 'v':
		case 'V':
			{
				printf("\n");
				if (theParam == MAX_SUBDIVS)
				{
					theParam = SCALE_FACT;
				}				
				else
				{
					theParam = MAX_SUBDIVS;
					printf("Scene  - Maximum Subdi(v)s: %i\r", g_scene->maxSubdivs());
					fflush(stdout);
				}
				break;
			}

		case 'n':
		case 'N':
			{
				printf("\n");
				if (theParam == NOISE_THRESHOLD)
				{
					theParam = SCALE_FACT;
				}				
				else
				{
					theParam = NOISE_THRESHOLD;
					printf("Scene  - (N)oise Threshold: %.6f\r", g_scene->noise());
					fflush(stdout);
				}
				break;
			}

        default:
			break;
    }
    glutPostRedisplay();
}


void
MiroWindow::reshape(int w, int h)
{
    g_image->resize(w, h);
    glViewport(0, 0, w, h);
    g_camera->setRenderer(Camera::RENDER_OPENGL);
    glutPostRedisplay();
}

