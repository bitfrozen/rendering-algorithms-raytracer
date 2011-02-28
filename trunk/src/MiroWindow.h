#ifndef CSE168_MIRO_GLUT_WINDOW_H_INCLUDED
#define CSE168_MIRO_GLUT_WINDOW_H_INCLUDED

class MiroWindow
{
public:
    MiroWindow(int * argc, char* argv[]);
    
    void mainLoop();
    
    // GLUT Event handlers
    void display();
    void reshape(int x, int y);
    void keyboard(unsigned char key, int x, int y);
	void specialKeys(int key, int x, int y);
    void mouse(int btn, int state, int x, int y);
    void motion(int x, int y);
	void increaseParam();
	void decreaseParam();
    
protected:
	enum paramType_t {SCALE_FACT, FOV, FOCUS_PLANE, APERTURE, NUM_PATHS, NUM_BOUNCES, MIN_SUBDIVS, MAX_SUBDIVS, NOISE_THRESHOLD, SHUTTER_SPEED};

	paramType_t theParam;

    float m_scaleFact;
	float m_microScaleFact;
    int   m_activeButton;
    int   m_mouseX, m_mouseY;
};

#endif // CSE168_MIRO_GLUT_WINDOW_H_INCLUDED


