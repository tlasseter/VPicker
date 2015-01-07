package com.Sts.Utilities;

import com.sun.opengl.util.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.nio.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsJOGLShader
{
    static public boolean initialized = false;
    static public boolean canUseShader = false;
    static public int nBytesPerTexel = 1;
    static public final int ARB_TLUT_NO_SPECULAR_LIGHTS = 1;
	static public final int ARB_TLUT_WITH_SPECULAR_LIGHTS = 2;
	static public final int ARB_TLUT_NO_SPECULAR_LIGHTS_ALPHA = 3;
	static public final int NONE = 0;
//    static public final int UNINITIALIZED = -1;

    static int RESERVED_TEXTURE_ID = 0;
	static private int activeType = 0;

    static String TLUT_2D_Fragment_Program =
        "!!ARBfp1.0\n" +
        "TEMP tTexCoord, tColor, t2Color,tBias;\n" +
        "PARAM bias = { -.0019 , -0.0019, -0.0019, -0.0019 };\n" +
        // Look up the texture coordinate from the fragment
        "TEX  tTexCoord, fragment.texcoord[0], texture[0], 2D;\n" +
        // bias a bit
        "ADD  tBias, tTexCoord.argb, bias;\n" +
        // dependent-texture read, use the texture coordinate as // index into the color table (2 1-D color tables stitched into a // 2D texture)
        "TEX  tColor, tBias.argb, texture[1], 2D;\n" +
        "MUL  result.color, tColor, fragment.color;\n" +
        "END\n";

	static String TLUT_2D_Fragment_Program_Alpha =
		"!!ARBfp1.0\n" +
		"TEMP tTexCoord, tColor, t2Color,tBias;\n" +
		"PARAM bias = { -.0019 , -0.0019, -0.0019, -0.0019 };\n" +
		// Look up the texture coordinate from the fragment
		"TEX  tTexCoord, fragment.texcoord[0], texture[0], 2D;\n" +
		// bias a bit
		"ADD  tBias, tTexCoord.argb, bias;\n" +
		// dependent-texture read, use the texture coordinate as // index into the color table (2 1-D color tables stitched into a // 2D texture)
		"TEX  tColor, tBias.argb, texture[1], 2D;\n" +
		"MUL  t2Color, tColor, fragment.color;\n" +
		// "MOV  result.color, t2Color;\n" +
		// tTexCoord.a is the alpha in a luminance-alpha data type
		"MUL  result.color, t2Color, tTexCoord.a;\n" +
        "MOV  result.color.a, t2Color.a;\n" +
        "END\n";

    static String TLUT_2D_Fragment_Program_Old =
		"!!ARBfp1.0\n" +
		"TEMP tTexCoord, tColor, t2Color;\n" +
		"TEX  tTexCoord, fragment.texcoord[0], texture[0], 2D;\n" +
		"TEX  tColor, tTexCoord.argb, texture[1], 2D;\n" +
		"MUL  result.color, tColor, fragment.color;\n" +
		"END\n";
	static String TLUT_2D_Light_Fragment_Program =
		"!!ARBfp1.0\n" +
		"TEMP tTexCoord, tColor, t2Color;\n" +
		"TEX  tTexCoord, fragment.texcoord[0], texture[0], 2D;\n" +
		"TEX  tColor, tTexCoord.argb, texture[1], 2D;\n" +
		"MUL  t2Color, tColor, fragment.color;\n" +
		"ADD  result.color, t2Color, fragment.color.secondary;\n" +
		"MOV  result.color.a, t2Color.a;\n" +
		"END\n";

    public static void initialize(GL gl)
    {
        if(initialized) return;
        canUseShader = StsJOGLShader.canDoARBShader(gl);
    }

    public static boolean canDoARBShader()
	{
        return canUseShader;
    }

    /** ARB shader needs these three conditions to be satisfied. */
    public static boolean canDoARBShader(GL gl)
	{
		if(!gl.isExtensionAvailable("GL_ARB_fragment_program")) return false;
        if(!gl.isExtensionAvailable("GL_ARB_multitexture")) return false;
	    return getMaxTextureUnits(gl) >= 2;
	}

    static private int getMaxTextureUnits(GL gl)
    {
        int[] buf = new int[1];
        gl.glGetIntegerv(GL.GL_MAX_TEXTURE_UNITS, buf, 0);
        return buf[0];
    }
/*
	public static boolean reloadTLUT(GL gl)
	{
		int error;

		// quietly don't use shader when picking
		int[] mode = new int[1];
		gl.glGetIntegerv(GL.GL_RENDER_MODE,mode,0);
		if (mode[0] == GL.GL_SELECT) return true;
		gl.glPixelTransferf(GL.GL_MAP_COLOR, 0);

        if(!checkInitializeARBColormap(gl)) return false;
        gl.glActiveTexture(GL.GL_TEXTURE1);
        gl.glBindTexture(GL.GL_TEXTURE_2D, RESERVED_TEXTURE_ID);
		if(colormapBuffer == null)
        {
            StsException.systemError(StsJOGLShader.class, "reloadTLUT", "colormapBuffer cannot be null.");
            return false;
        }
        colormapBuffer.rewind();
        gl.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA8, 256, 2, 0, GL.GL_RGBA, GL.GL_FLOAT, colormapBuffer); //colormap
        error = gl.glGetError();
        if(error != 0)
        {
            GLU glu = new GLU();
            System.out.println("shader TLUT reload err code " + error + " " + glu.gluErrorString(error));
        }
        // back to texture 0
        gl.glActiveTexture(GL.GL_TEXTURE0);
    }
*/
    static private boolean isPicking(GL gl)
    {
        int[] mode = new int[1];
		gl.glGetIntegerv(GL.GL_RENDER_MODE,mode,0);
		return mode[0] == GL.GL_SELECT;
    }

    static public boolean loadEnableARBColormap(GL gl, float[][] rgbaArray, int type)
    {
        if(rgbaArray == null) return false;
        int nColors = rgbaArray[0].length;
        FloatBuffer colormapBuffer = computeColormapBuffer(rgbaArray, nColors);
        return loadEnableARBColormap(gl, colormapBuffer, type, nColors);
    }

    static public boolean loadEnableARBColormap(GL gl, FloatBuffer colormapBuffer, int type)
    {
        return loadEnableARBColormap(gl, colormapBuffer, type, 256);
    }

    static public boolean loadEnableARBColormap(GL gl, FloatBuffer colormapBuffer, int type, int colormapLength)
    {
        if(!loadARBColormap(gl, colormapBuffer, colormapLength)) return false;
        if(!enableARBShader(gl, type)) return false;
        return true;
    }

    static public boolean loadARBColormap(GL gl, FloatBuffer colormapBuffer, int colormapLength)
    {
        int error;
		error = gl.glGetError(); // clean error status
		gl.glPixelTransferf(GL.GL_MAP_COLOR, 0);

        if(!checkInitializeARBColormap(gl)) return false;

        if(colormapBuffer == null)
        {
            StsException.systemError(StsJOGLShader.class, "loadARBColormap", "colormapBuffer cannot be null.");
            return false;
        }
        // colormap goes into texture unit #1
        colormapBuffer.rewind();
        gl.glActiveTexture(GL.GL_TEXTURE1);
		gl.glBindTexture(GL.GL_TEXTURE_2D, RESERVED_TEXTURE_ID);
		gl.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA8, colormapLength, 2, 0, GL.GL_RGBA, GL.GL_FLOAT, colormapBuffer); //colormap
        // back to texture 0
		gl.glActiveTexture(GL.GL_TEXTURE0);
		error = gl.glGetError();
		if(error != 0)
		{
			GLU glu = new GLU();
			System.out.println("shader TLUT load err code " + error + " " + glu.gluErrorString(error));
		}
		return true;
	}

    static public boolean checkInitializeARBColormap(GL gl)
	{
		if(!canUseShader)return true;
		int error;
		error = gl.glGetError(); // clean error status

		gl.glPixelTransferf(GL.GL_MAP_COLOR, 0);

		// colormap goes into texture unit #1
		gl.glActiveTexture(GL.GL_TEXTURE1);
		gl.glEnable(GL.GL_TEXTURE_2D);
		gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE);
		gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE);
		gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST);
		gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST);

		if(RESERVED_TEXTURE_ID <= 0)
		{
			int[] textures = new int[1];
			gl.glGenTextures(1, textures, 0);
			RESERVED_TEXTURE_ID = textures[0];
		}
		// back to texture 0
		gl.glActiveTexture(GL.GL_TEXTURE0);
		error = gl.glGetError();
		if(error != 0)
		{
			GLU glu = new GLU();
			System.out.println("shader TLUT load err code " + error + " " + glu.gluErrorString(error));
		}
		return true;
	}

    static public boolean disableARBShader(GL gl)
	{
		return enableARBShader(gl, NONE);
	}

	static public boolean enableARBShader(GL gl, int type)
	{
		if(!canUseShader)return false;
		int error;
		//GLU glu = new GLU();
		error = gl.glGetError();

	// quietly don't use shader when picking
    	int[] mode = new int[1];
        gl.glGetIntegerv(GL.GL_RENDER_MODE,mode,0);
        if (mode[0] == GL.GL_SELECT) return true;

		if(type == ARB_TLUT_NO_SPECULAR_LIGHTS)
		{
			activeType = type;
            nBytesPerTexel = 1;
            gl.glBindProgramARB(GL.GL_FRAGMENT_PROGRAM_ARB, ARB_TLUT_NO_SPECULAR_LIGHTS);

			gl.glProgramStringARB(GL.GL_FRAGMENT_PROGRAM_ARB, GL.GL_PROGRAM_FORMAT_ASCII_ARB,
								  TLUT_2D_Fragment_Program.length(), TLUT_2D_Fragment_Program);

			gl.glEnable(GL.GL_FRAGMENT_PROGRAM_ARB);

			gl.glActiveTexture(GL.GL_TEXTURE1);
			gl.glEnable(GL.GL_TEXTURE_2D);
			gl.glActiveTexture(GL.GL_TEXTURE0);
			//gl.glEnable(GL.GL_TEXTURE_2D);
		}
		else if(type == ARB_TLUT_WITH_SPECULAR_LIGHTS)
		{
			activeType = type;
            nBytesPerTexel = 1;
            gl.glBindProgramARB(GL.GL_FRAGMENT_PROGRAM_ARB, ARB_TLUT_WITH_SPECULAR_LIGHTS);

			gl.glProgramStringARB(GL.GL_FRAGMENT_PROGRAM_ARB, GL.GL_PROGRAM_FORMAT_ASCII_ARB,
								  TLUT_2D_Light_Fragment_Program.length(), TLUT_2D_Light_Fragment_Program);

			gl.glEnable(GL.GL_FRAGMENT_PROGRAM_ARB);

			//gl.glActiveTexture(GL.GL_TEXTURE1);
			//gl.glEnable(GL.GL_TEXTURE_2D);
			gl.glActiveTexture(GL.GL_TEXTURE0);
			//gl.glBindTexture(GL.GL_TEXTURE_2D, 0);
			//gl.glEnable(GL.GL_TEXTURE_2D);
		}
		else if(type == ARB_TLUT_NO_SPECULAR_LIGHTS_ALPHA)
		{
			activeType = type;
            nBytesPerTexel = 2;
            gl.glBindProgramARB(GL.GL_FRAGMENT_PROGRAM_ARB, ARB_TLUT_NO_SPECULAR_LIGHTS_ALPHA);
			gl.glProgramStringARB(GL.GL_FRAGMENT_PROGRAM_ARB, GL.GL_PROGRAM_FORMAT_ASCII_ARB,
								  TLUT_2D_Fragment_Program_Alpha.length(), TLUT_2D_Fragment_Program_Alpha);

			gl.glEnable(GL.GL_FRAGMENT_PROGRAM_ARB);

			gl.glActiveTexture(GL.GL_TEXTURE1);
			gl.glEnable(GL.GL_TEXTURE_2D);
			gl.glActiveTexture(GL.GL_TEXTURE0);
			//gl.glEnable(GL.GL_TEXTURE_2D);
		}
		else
		{
			gl.glDisable(GL.GL_FRAGMENT_PROGRAM_ARB);
			gl.glActiveTexture(GL.GL_TEXTURE1);
			gl.glDisable(GL.GL_TEXTURE_2D);
			gl.glActiveTexture(GL.GL_TEXTURE2);
			gl.glDisable(GL.GL_TEXTURE_2D);
			gl.glActiveTexture(GL.GL_TEXTURE0);
			//gl.glDisable(GL.GL_TEXTURE_2D);
			int[] buf = new int[1];
			buf[0] = activeType;
			gl.glDeleteProgramsARB(1, buf, 0);
		}

		error = gl.glGetError();
		if(error != 0)
		{
			GLU glu = new GLU();
			System.out.println("shader " + type + " err code " + error + " " + glu.gluErrorString(error));
		}
		return true;
	}

    static public void close()
    {
        initialized = false;
        RESERVED_TEXTURE_ID = 0;
    }

    static public FloatBuffer computeColormapBuffer(float[][] arrayRGBA, int nColors)
    {
        FloatBuffer colormapBuffer = BufferUtil.newFloatBuffer(nColors * 4 * 2);
		colormapBuffer.rewind();
		for(int pass = 0; pass < 2; pass++)
        {
            for(int i = 0; i < nColors; i++)
			{
				for(int j = 0; j < 4; j++)
				{
					colormapBuffer.put(arrayRGBA[j][i]);
				}
			}
        }
        return colormapBuffer;
    }
}
