package com.Sts.Utilities.Shaders;

import com.Sts.Utilities.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 17, 2009
 * Time: 6:21:21 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsUniformVariableF extends StsUniformVariable
{
    float[] values;

    StsUniformVariableF(String name, float[] values)
    {
        super(name, values.length);
        this.values = values;
    }

    StsUniformVariableF(String name, float value)
    {
        super(name, 1);
        this.values = new float[] { value };
    }

    void setValues(GL gl, int handle, float[] values, int length)
    {
        if(length == dimension)
            this.values = values;
        else
        {
            length = Math.min(dimension, length);
            System.arraycopy(values, 0, this.values, 0, length);
        }
        setValues(gl, handle);
    }

    void setValues(GL gl, int handle)
    {
        if(location < 0)
        {
            location = getLocation(gl, handle);
		    if (location < 0)
                System.out.println("bad location for "+name);
        }
        switch(dimension)
        {
            case 1:
                gl.glUniform1f(location, values[0]);
                break;
            case 2:
                gl.glUniform2f(location, values[0], values[1]);
                break;
            case 3:
                gl.glUniform3f(location, values[0], values[1], values[2]);
                break;
            case 4:
                gl.glUniform4f(location, values[0], values[1], values[2], values[3]);
                break;
            default:
                StsException.systemError(StsBeachballShader.class, "setValue", "Called with undefined dimension: " + dimension);
        }
    }
}
