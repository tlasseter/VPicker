package com.Sts.Utilities;

import javax.media.opengl.*;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author Tom Lasseter
 * @version 1.0
 */

public interface StsGLPicker
{
    public boolean execute(GLEventListener caller, GLAutoDrawable glDrawable) throws Exception;

    public boolean doPick(GLAutoDrawable glDrawable);
}
