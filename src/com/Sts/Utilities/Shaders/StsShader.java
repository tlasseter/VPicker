package com.Sts.Utilities.Shaders;

import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.lang.reflect.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 14, 2009
 * Time: 10:54:32 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsShader
{
    Class shaderClass;
    String shortClassname;
    int vertexShader = -1;
    int fragmentShader = -1;
    int shaderProgram = -1;

    static HashMap<Class, StsShader> shaders = new HashMap();
    static StsShader currentShader = null;

    static int[] operationOk = new int[1];

    public StsShader()
    {
    }

    public StsShader(GL gl, String shaderName) throws GLException
    {
        if(!gl.glGetString(GL.GL_EXTENSIONS).contains("GL_ARB_shading_language_100"))
        {
            systemError(gl, null, "getEnableShader", "GL_ARB shader language not available.");
            throw new GLException("GL_ARB shader language not available");
        }

        this.shaderClass = getClass();
        shortClassname = StsToolkit.getSimpleClassname(shaderClass);
        vertexShader = gl.glCreateShaderObjectARB(GL.GL_VERTEX_SHADER);
        fragmentShader = gl.glCreateShaderObjectARB(GL.GL_FRAGMENT_SHADER);
        String vertexShaderString = getShaderSource(shaderName + "VertexShader");
        gl.glShaderSourceARB(vertexShader, 1, new String[] { vertexShaderString}, null, 0);
        String fragmentShaderString = getShaderSource(shaderName + "FragmentShader");
        gl.glShaderSourceARB(fragmentShader, 1, new String[] { fragmentShaderString}, null, 0);
        gl.glCompileShaderARB(vertexShader);
        gl.glGetObjectParameterivARB(vertexShader, GL.GL_OBJECT_COMPILE_STATUS_ARB, operationOk, 0);
        if(operationOk[0] != GL.GL_TRUE)
        {
            systemError(gl, this, shaderName + " constructor", "Failed to compile vertex shader");
            return;
        }
        gl.glCompileShaderARB(fragmentShader);
        gl.glGetObjectParameterivARB(fragmentShader, GL.GL_OBJECT_COMPILE_STATUS_ARB, operationOk, 0);
        if(operationOk[0] != GL.GL_TRUE)
        {
            systemError(gl, this, shaderName + " constructor", "Failed to compile fragment shader");
            return;
        }

        shaderProgram = gl.glCreateProgramObjectARB();
        gl.glAttachObjectARB(shaderProgram, vertexShader);
        gl.glDeleteObjectARB(vertexShader);
        vertexShader = -1;
        gl.glAttachObjectARB(shaderProgram, fragmentShader);
        gl.glDeleteObjectARB(fragmentShader);
        fragmentShader = -1;
        gl.glLinkProgramARB(shaderProgram);
        gl.glValidateProgramARB(shaderProgram);
        StsException.systemDebug(StsShader.class, shaderName + " constructor", "info log string: " + getInfoLog(gl));
        gl.glGetObjectParameterivARB(shaderProgram, GL.GL_OBJECT_LINK_STATUS_ARB, operationOk, 0);
        if(operationOk[0] != GL.GL_TRUE)
        {
            systemError(gl, this, "constructor", "Failed to link " + shaderName + " shader program");
            return;
        }
    }

    String getSimpleClassname() { return StsToolkit.getSimpleClassname(shaderClass); }

    static StsShader getEnableShader(GL gl, Class shaderClass)
    {
        if(currentShader != null && currentShader.shaderClass == shaderClass) return currentShader;
        StsShader shader = shaders.get(shaderClass);
        if(shader == null)
        {
            try
            {
                Constructor constructor = shaderClass.getConstructor(GL.class );
                shader = (StsShader)constructor.newInstance(gl);
                shaders.put(shaderClass, shader);

            }
            catch(Exception e)
            {
                StsException.outputWarningException(StsShader.class, "getEnableShader", "Failed to find shader " + StsToolkit.getSimpleClassname(shaderClass), e);
                return null;
            }
        }
        currentShader = shader;
        gl.glUseProgramObjectARB(shader.shaderProgram);
        return shader;
    }

    static public void disableCurrentShader(GL gl)
    {
        if(currentShader == null) return;
        currentShader.disableShader(gl);
        currentShader = null;
    }

    private void disableShader(GL gl)
     {
         gl.glUseProgramObjectARB(0);
     }

    static public void deleteProgram(Class shaderClass)
    {
        shaders.remove(shaderClass);
    }

    static void systemError(GL gl, StsShader shader, String methodName, String message)
    {
        String info = "no info";
        if(shader != null) info = shader.getInfoLog(gl);
        StsException.systemError(StsShader.class, shader.shortClassname, "method: " + methodName + " " + message + " info: " + info);
    }

    static String getShaderSource(String shaderName)
    {
        String source = "";
        BufferedReader bufferedReader = null;
        try
        {
            java.net.URL url = StsShader.class.getResource("resources/" + shaderName + ".txt");
            InputStream inputStream = url.openStream();
            bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
        }
        catch(Exception e)
        {
            StsException.systemError(StsShader.class, "getShaderSource", "Failed to construct inputstream from shader source: " + shaderName);
            return source;
        }

        String line;
        try
        {
            while ((line=bufferedReader.readLine()) != null)
                source += line + "\n";
            return source;
        }
        catch(Exception e)
        {
            StsException.systemError(StsShader.class, "readShaderSource", "Failed reading shader source: " + shaderName);
            return source;
        }
    }

    public String getInfoLog(GL gl)
    {
        int[] length = new int[1];
        gl.glGetObjectParameterivARB(shaderProgram, GL.GL_OBJECT_INFO_LOG_LENGTH_ARB, length, 0);
        int lengthBytes = length[0];
        if(lengthBytes == 0) return "";
        byte[] infoLog = new byte[lengthBytes];
        gl.glGetInfoLogARB(shaderProgram, 2056, length, 0, infoLog, 0);
        gl.glGetObjectParameterivARB(shaderProgram, GL.GL_OBJECT_VALIDATE_STATUS_ARB, operationOk, 0);
        if(operationOk[0] != GL.GL_TRUE)
            return "Failed to validate.";
        return new String(infoLog);
    }
}
