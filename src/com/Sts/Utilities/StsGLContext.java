
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;

public interface StsGLContext
{
  GL getGL();
  GLU getGLU();
}

