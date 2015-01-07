package com.Sts.Utilities.Shaders;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 17, 2009
 * Time: 6:20:10 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsUniformVariable
{
    int location = -1;
    String name;
    int dimension;

    abstract void setValues(GL gl, int handle);

    StsUniformVariable()
    {
    }

    StsUniformVariable(String name, int dimension)
    {
        this.name = name;
        this.dimension = dimension;
    }

    int getLocation(GL gl, int handle)
    {
        return gl.glGetUniformLocationARB(handle, name);
    }
}
