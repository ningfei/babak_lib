#ifndef _landmarks_h

float detect_lm(SPH &searchsph, SPH &testsph, SHORTIM testim, int lmcm[], SPH &refsph, int lm[]);
float *detect_landmarks(const char *subfile, const char *mdlfile, int &nl, char ppmflg);
float *detect_landmarks(SHORTIM subim, const char *mdlfile, int &nl);
float *read_landmark_centers(const char *mdlfile, int &nl);

#define _landmarks_h

#endif
