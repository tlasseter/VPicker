
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.io.*;
import java.util.*;

/**
  * Color is defined by 4 floats for rgba in the range 0.0f to 1.0f
  * @author TJLasseter
  */

public class StsColor extends StsSerialize implements Cloneable, Serializable
{
    /** the rgb values as normalized floats (0.0f to 1.0f) */
	public final float red, green, blue, alpha;
    /** index of this color in a spectrum or vector of colors */
    public int idx = -1;
    // Colors in Comboboxes
    static public final StsColor RED = new StsColor(1.0f, 0.0f, 0.0f);
    static public final StsColor GREEN = new StsColor(0.0f, 1.0f, 0.0f);
    static public final StsColor BLUE = new StsColor(0.0f, 0.0f, 1.0f);
    static public final StsColor YELLOW = new StsColor(1.0f, 1.0f, 0.0f);
    static public final StsColor DARKKHAKI = new StsColor(189,183,107);
    static public final StsColor MAGENTA = new StsColor(1.0f, 0.0f, 1.0f);
    static public final StsColor CYAN = new StsColor(0.0f, 1.0f, 1.0f);
    static public final StsColor VIOLET = new StsColor(79,47,79);
    static public final StsColor CORAL = new StsColor(238,106,80);
    static public final StsColor DARKOLIVE = new StsColor(85,107,47);
    static public final StsColor BRONZE = new StsColor(166,125,61);
    static public final StsColor BROWN = new StsColor(166,42,42);
    static public final StsColor PURPLE = new StsColor(128,0,128);
    static public final StsColor ORANGE = new StsColor(255,165,0);
    static public final StsColor GRAY = new StsColor(0.50f, 0.50f, 0.50f);
    static public final StsColor BRASS = new StsColor(181,166,66);
    static public final StsColor DARK_RED = new StsColor(0.65f, 0.0f, 0.0f);
    static public final StsColor DARK_GREEN = new StsColor(0.0f, 0.65f, 0.0f);
    static public final StsColor DARK_BLUE= new StsColor(0.0f, 0.0f, 0.65f);
    static public final StsColor DARK_YELLOW = new StsColor(0.65f, 0.65f, 0.0f);
    static public final StsColor PINK = new StsColor(0.98f, 0.04f, 0.7f);
    static public final StsColor DARK_MAGENTA = new StsColor(0.65f, 0.0f, 0.65f);
    static public final StsColor MAROON = new StsColor(128,0,0);
    static public final StsColor DARK_CYAN = new StsColor(0.0f, 0.65f, 0.65f);
    static public final StsColor FUCHSIA = new StsColor(255,0,255);
    static public final StsColor COPPER = new StsColor(184,115,51);
    static public final StsColor DARKPURPLE = new StsColor(135,31,120);
    static public final StsColor DARK_GRAY = new StsColor(0.25f, 0.25f, 0.25f);
    static public final StsColor YELLOWGREEN = new StsColor(153,204,50);
    static public final StsColor BLACK = new StsColor(0.0f, 0.0f, 0.0f);
    static public final StsColor WHITE = new StsColor(1.0f, 1.0f, 1.0f);
    static public final StsColor GREY = new StsColor(0.50f, 0.50f, 0.50f);

    static public final StsColor DARK_GREY = new StsColor(0.25f, 0.25f, 0.25f);
    static public final StsColor LIGHT_GREY = new StsColor(0.75f, 0.75f, 0.75f);
    static public final StsColor LIGHT_GRAY = new StsColor(0.75f, 0.75f, 0.75f);
    static public final StsColor CLEAR = new StsColor(0.0f, 0.0f, 0.0f, 0.0f);
    static public final StsColor ALICEBLUE = new StsColor(240, 248, 255);
    static public final StsColor BLUEVIOLET = new StsColor(138, 43, 226);
    static public final StsColor CADETBLUE = new StsColor(95, 159, 159);
    static public final StsColor DARKTURQUOISE = new StsColor(0,206,209);
    static public final StsColor DARKSLATEBLUE = new StsColor(72,61,139);
    static public final StsColor DEEPSKYBLUE = new StsColor(0,191,255);
    static public final StsColor LIGHTBLUE = new StsColor(173,216,230);
    static public final StsColor LIGHTSLATEBLUE = new StsColor(132,112,255);
    static public final StsColor AQUAMARINE= new StsColor(112,219,147);
    static public final StsColor MIDNIGHTBLUE = new StsColor(25,25,112);
    static public final StsColor ROYALBLUE = new StsColor(65,105,225);
    static public final StsColor SLATEBLUE = new StsColor(106,90,205);
    static public final StsColor AZURE = new StsColor(193,205,205);
    // Browns
    static public final StsColor TEAL = new StsColor(0,128,128);
    static public final StsColor ROSYBROWN = new StsColor(188,143,143);
    static public final StsColor SADDLEBROWN = new StsColor(139,69,19);
    static public final StsColor SANDYBROWN = new StsColor(244,164,96);
    static public final StsColor DARKBROWN = new StsColor(92,64,51);
    static public final StsColor BURLYWOOD = new StsColor(222,184,135);
    static public final StsColor TAN = new StsColor(210,180,140);
    static public final StsColor DARKTAN = new StsColor(151,105,79);
    static public final StsColor SIENNA = new StsColor(142,107,35);
    // Greens
    static public final StsColor OLIVE = new StsColor(128,128,0);
    static public final StsColor LIGHTOLIVE = new StsColor(162,205,90);
    static public final StsColor DARKSEAGREEN = new StsColor(143,188,143);
    static public final StsColor FORESTGREEN = new StsColor(34,139,34);
    static public final StsColor LAWNGREEN = new StsColor(124,252,0);
    static public final StsColor LIMEGREEN = new StsColor(50,205,50);
    static public final StsColor SPRINGGREEN = new StsColor(0,250,154);
    static public final StsColor PALEGREEN = new StsColor(142,107,35);
    static public final StsColor SEAGREEN = new StsColor(84,255,159);
    static public final StsColor CHARTREUSE = new StsColor(118,238,0);
    //Oranges
    static public final StsColor DARKORANGE = new StsColor(238,118,0);
    static public final StsColor DARKSALMON = new StsColor(255,160,122);
    static public final StsColor LIGHTCORAL = new StsColor(240,128,128);
    static public final StsColor PEACHPUFF = new StsColor(255,218,185);
    static public final StsColor BISQUE = new StsColor(238,213,183);
    static public final StsColor MANDARIANORANGE = new StsColor(142,35,35);
    static public final StsColor ORANGERED = new StsColor(255,36,0);
    // Reds
    static public final StsColor DEEPPINK = new StsColor(255,20,147);
    static public final StsColor HOTPINK = new StsColor(255,105,180);
    static public final StsColor INDIANRED = new StsColor(205,85,85);
    static public final StsColor MEDIUMVIOLETRED = new StsColor(199,21,133);
    static public final StsColor PALEVIOLETRED = new StsColor(219,112,147);
    static public final StsColor VILOETRED = new StsColor(208,32,144);
    static public final StsColor FIREBRICK = new StsColor(205,38,38);
    static public final StsColor FLESH = new StsColor(245,204,176);
    static public final StsColor FELDSPAR = new StsColor(209,146,117);
    static public final StsColor TOMATO = new StsColor(255,99,71);
    static public final StsColor DUSTYROSE = new StsColor(133,99,99);
    static public final StsColor SPICYPINK = new StsColor(255,28,174);
    // Violets
    static public final StsColor PLUM = new StsColor(234,173,234);
    static public final StsColor VIOLETBLUE = new StsColor(159,95,159);
    static public final StsColor THISTLE = new StsColor(205,181,205);
    static public final StsColor ORCHID = new StsColor(219,112,219);
    static public final StsColor MEDIUMPURPLE = new StsColor(159,121,238);
    static public final StsColor DARKORCHID = new StsColor(173,58,238);
    // Whites
    static public final StsColor ANTIQUEWHITE = new StsColor(238,223,204);
    static public final StsColor GHOSTWHITE = new StsColor(248,248,255);
    static public final StsColor NAVAJOWHITE = new StsColor(255,222,173);
    static public final StsColor IVORY = new StsColor(255,255,240);
    static public final StsColor SEASHELL = new StsColor(205,197,191);
    static public final StsColor WHEAT = new StsColor(238,216,174);
    // Yellows
    static public final StsColor DARKGOLDENROD = new StsColor(255,185,15);
    static public final StsColor LEMONCHIFFON = new StsColor(255,250,205);
    static public final StsColor LIGHTGOLDENROD = new StsColor(255,236,139);
    static public final StsColor LIGHTYELLOW = new StsColor(238,232,170);
    static public final StsColor GOLDENROD = new StsColor(218,165,32);
    static public final StsColor GOLDEN = new StsColor(238,201,0);
    // Metallics
    static public final StsColor GOLD = new StsColor(205,127,50);
    static public final StsColor SILVER = new StsColor(192,192,192);
    static public final StsColor STEELBLUE = new StsColor(35,107,142);

    static public final String[] colorNames32 = {"Red","Green","Blue","Yellow","Khaki","Magenta","Cyan","Violet",
                          "Coral","Olive", "Bronze", "Brown","Purple","Orange","Gray","Brass","Dark Red", "Dark Green", "Dark Blue",
                           "Dark Yellow","Pink", "Dark Magenta", "Maroon", "Dark Cyan", "Fuchsia", "Copper", "Dark Purple",
                           "Dark Grey", "Yellow Green", "Black", "White", "Grey"};
    static public final StsColor[] colors32 = {RED, GREEN, BLUE, YELLOW, DARKKHAKI, MAGENTA, CYAN, VIOLET, CORAL,
            DARKOLIVE, BRONZE, BROWN, PURPLE, ORANGE, GRAY, BRASS, DARK_RED, DARK_GREEN, DARK_BLUE, DARK_YELLOW, PINK, DARK_MAGENTA,
            MAROON, DARK_CYAN, FUCHSIA, COPPER, DARKPURPLE, DARK_GRAY, YELLOWGREEN, BLACK, WHITE, GREY };

    static public final String[] bkgdColorNames = {"BLACK","DARK_GREY","LIGHT_GREY","WHITE","KHAKI","COPPER","BRONZE","BRASS","SILVER"};
    static public final StsColor[] bkgdColors = {BLACK,DARK_GREY,LIGHT_GREY,WHITE,DARKKHAKI,COPPER,BRONZE,BRASS,SILVER};

    static public final String[] colorNames = {"BLACK","RED","GREEN","BLUE","YELLOW","MAGENTA","CYAN","WHITE","GREY",
                                               "DARK_GREY","LIGHT_GREY","GRAY","DARK_GRAY","LIGHT_GRAY","ORANGE","PINK","PURPLE","BROWN", "DARK_GREEN"};
    static public final StsColor[] colors = {BLACK,RED,GREEN,BLUE,YELLOW,MAGENTA,CYAN,WHITE,GREY,
                                             DARK_GREY,LIGHT_GREY,GRAY,DARK_GRAY,LIGHT_GRAY,ORANGE,PINK,PURPLE,BROWN, DARK_GREEN};

    static public final StsColor[] colors8 = new StsColor[] { RED, GREEN, BLUE, YELLOW, MAGENTA, CYAN,  ORANGE, PINK };
    static public final String[] colorNames8 = new String[] {"RED", "GREEN", "BLUE", "YELLOW", "MAGENTA", "CYAN", "ORANGE", "PINK" };

    static public final StsColor[] colors16 = new StsColor[] { RED, GREEN, BLUE, YELLOW, MAGENTA, CYAN,  ORANGE, PINK, PURPLE, BROWN, DARK_RED, DARK_GREEN, DARK_BLUE, DARK_YELLOW, DARK_MAGENTA, DARK_CYAN  };
    static public final String[] colorNames16 = new String[] {"RED", "GREEN", "BLUE", "YELLOW", "MAGENTA", "CYAN", "ORANGE", "PINK", "PURPLE", "BROWN", "DARK_RED", "DARK_GREEN", "DARK_BLUE", "DARK_YELLOW", "DARK_MAGENTA", "DARK_CYAN" };

    static final float norm = 1.0f/255.0f;

    /** scratch vector used for 4-component color */
    static float[] rgba = new float[4];

    /** constructor for db operations */
    public StsColor()
    {
        red = 0.0f; green = 0.0f; blue = 0.0f; alpha = 1.0f;
    }

   	/** Constructs and initializes a StsColor from the rgb normalized values. */
   	public StsColor(float red, float green, float blue)
	{
         this.red = red; this.green = green; this.blue = blue; alpha = 1.0f;
	}

    public StsColor(float red, float green, float blue, float alpha)
    {
        this.red = red; this.green = green; this.blue = blue; this.alpha = alpha;
	}

    public StsColor(int idx, float red, float green, float blue, float alpha)
    {
        this.idx = idx;
        this.red = red; this.green = green; this.blue = blue; this.alpha = alpha;
	}

   	/** Constructs and initializes a StsColor from the rgb normalized vector. */
       public StsColor(float[] rgb)
	{
         if(rgb.length == 4)
         {
             this.red = rgb[0];
             this.green = rgb[1];
             this.blue = rgb[2];
             this.alpha = rgb[3];
         }
         else if(rgb.length == 5)
         {
             this.idx = (int)rgb[0];
             this.red = rgb[1];
             this.green = rgb[2];
             this.blue = rgb[3];
             this.alpha = rgb[4];
         }
         else
         {
             this.red = rgb[0];
             this.green = rgb[1];
             this.blue = rgb[2];
             this.alpha = 1.0f;
         }
	}

   	/** Constructs and initializes a StsColor from a standard java Color. */
	public StsColor(Color color)
	{
//        this(norm*color.getRed(), norm*color.getGreen(), norm*color.getBlue(), norm*color.getAlpha());
        float[] rgba = new float[4];
        color.getComponents(rgba);
        red = rgba[0];
        green = rgba[1];
        blue = rgba[2];
        alpha = rgba[3];
    }

    /** Constructs and initializes a StsColor from a standard java Color and index. */
    public StsColor(Color color, int idx)
    {
        float[] rgba = new float[4];
        color.getComponents(rgba);
        red = rgba[0];
        green = rgba[1];
        blue = rgba[2];
        alpha = rgba[3];
        this.idx = idx;
    }

   	/** Colors are 0 to 255; normalize to 0.0 to 1.0. */
	public StsColor(int red, int green, int blue)
	{
        this(norm*red, norm*green, norm*blue);
    }

    public StsColor(int red, int green, int blue, int alpha)
    {
        this(norm*red, norm*green, norm*blue, norm*alpha);
    }

    public StsColor(int idx, int red, int green, int blue, int alpha)
    {
        this(idx, norm*red, norm*green, norm*blue, norm*alpha);
    }

    public boolean initialize(StsModel model) { return true; }

  	/** Constructs and initializes a StsColor from an StsColor. */
	public StsColor(StsColor color)
	{
        this(color.idx, color.red, color.green, color.blue, color.alpha);
    }

    public StsColor(StsColor color, float alpha)
    {
         this(color.idx, color.red, color.green, color.blue, alpha);
    }

    final public void setGLColor(GL gl)
    {
        gl.glColor4f(red, green, blue, alpha);
    }

    static public final void setGLColor(GL gl, float red, float green, float blue, float alpha)
    {
        gl.glColor4f(red, green, blue, alpha);
    }

    static public final void setGLJavaColor(GL gl, Color color)
	{
        color.getComponents(rgba);
		gl.glColor4fv(rgba, 0);
	}

   /**
     * Makes an StsColor from normalized hue, saturation, brightness components.
     * @param h hue(0.0f to 1.0f)
     * @param s saturation (0.0f to 1.0f)
     * @param b brightness (0.0f to 1.0f)
     * @return the StsColor
     */
    final public StsColor getHSBColor(float h, float s, float b) throws StsException
    {
    	Color color = Color.getHSBColor(h, s, b);
        float red = norm*color.getRed();
        float green = norm*color.getGreen();
        float blue = norm*color.getBlue();
        return new StsColor(red, green, blue);
    }

    /** Returns a Color from this StsColor
      * @return the Color
      */
    final public Color getColor()
    {
    	return new Color(red, green, blue, alpha);
    }
	/** Returns a Color from this StsColor
	  * @return the Color
	  */
	static final public Color getInverseColor(StsColor c)
	{
		return new Color(1.f-c.red, 1.f-c.green, 1.f-c.blue, c.alpha);
	}
	static final public StsColor getInverseStsColor(StsColor c)
		{
			return new StsColor(1.f-c.red, 1.f-c.green, 1.f-c.blue, c.alpha);
		}


    static public final Color[] getColors(StsColor[] stsColors)
    {
        int nColors = stsColors.length;
        Color[] colors = new Color[nColors];
        for(int n = 0; n < nColors; n++)
            colors[n] = stsColors[n].getColor();
        return colors;
    }

    final public String getHSBString()
    {
        float[] hsbvals = new float[3];
        Color color = new Color(red, green, blue);
        color.RGBtoHSB(color.getRed(), color.getGreen(), color.getBlue(), hsbvals);

        int hue = (int)(360f*hsbvals[0]);
        int bright = (int)(100f*hsbvals[1]);
        int sat = (int)(100f*hsbvals[2]);

        return new String("h: " + hue + " s: " + (float)bright/100f + " b: " + (float)sat/100f);
    }
/*
    final public void setBeachballColors(Color color)
    {
        if(color == null) return;
        float[] rgba = new float[4];
        color.getComponents(rgba);
        setBeachballColors(rgba);
    }

    final public void setBeachballColors(float[] rgba)
    {
        red = rgba[0];
        green = rgba[1];
        blue = rgba[2];
        alpha = rgba[3];
//        currentModel.instanceChange(this, "colorChanged");
    }

    final public void setBeachballColors(int red, int green, int blue, int alpha)
    {
        this.red = norm*red;
        this.green = norm*green;
        this.blue = norm*blue;
        this.alpha = norm*alpha;
//        currentModel.instanceChange(this, "colorChanged");
    }
*/
    /** get RGB values in array */
    final public float[] getRGB()
    {
        float[] rgb = new float[3];
        rgb[0] = red;
        rgb[1] = green;
        rgb[2] = blue;
        return rgb;
    }

    /** get RGBA values in array */
    final public float[] getRGBA()
    {
        return new float[] { red, green, blue, alpha };
    }

    public float[][] getRGBAArray()
	{
        return new float[][] { { red }, { green }, { blue }, { alpha } };
    }
    /** get iRGBA values in array */
    final public float[] getIRGBA()
    {
        float[] irgba = new float[5];
        irgba[0] = (float)idx;
        irgba[1] = red;
        irgba[2] = green;
        irgba[3] = blue;
        irgba[4] = alpha;
        return irgba;
    }

    public boolean rgbEquals(Color color)
    {
        if(Math.round(red*255) != color.getRed()) return false;
        if(Math.round(green*255) != color.getGreen()) return false;
        if(Math.round(blue*255) != color.getBlue()) return false;
        return true;
    }
    static public StsColor getColorByName32(String name)
    {
        for(int i=0; i< colorNames32.length; i++)
            if(colorNames32[i].equalsIgnoreCase(name))
                return colors32[i];
        return RED;
    }
    static public String getNameByStsColor32(StsColor color)
    {
        for(int i=0; i< colorNames32.length; i++)
            if(colors32[i].equals(color))
                return colorNames32[i];
        return "RED";
    }
    static public StsColor getColorByName16(String name)
    {
        for(int i=0; i< colorNames16.length; i++)
            if(colorNames16[i].equalsIgnoreCase(name))
                return colors16[i];
        return RED;
    }
    static public String getNameByStsColor16(StsColor color)
    {
        for(int i=0; i< colorNames16.length; i++)
            if(colors16[i].equals(color))
                return colorNames16[i];
        return "RED";
    }
    static public StsColor getColorByName(String name)
    {
        for(int i=0; i<colorNames.length; i++)
            if(colorNames[i].equalsIgnoreCase(name))
                return colors[i];
        return RED;
    }
    static public String getNameByColor(Color color)
    {
        for(int i=0; i<colors.length; i++)
            if(colors[i].getColor().equals(color))
                return colorNames[i];
        return "RED";
    }

    public String toLabelString()
    {
        return new String(red + " " + green + " " + blue + " " + alpha);
    }

    public String toString()
	{
		return "r " + red + " g " + green + " b " + blue + " a " + alpha;
	}

    static public int getColorIndex(StsColor color, StsColor[] colors)
    {
        for(int n = 0; n < colors.length; n++)
            if(colors[n] == color) return n;
        return -1;
    }

    static public StsColor colorFromString(String string)
    {
        if(string == null) return StsColor.BLACK;
        StringTokenizer stok = new StringTokenizer(string, " ");
	    int nTokens = stok.countTokens();
        if(nTokens == 1)
        {
            Color color = Color.decode(string);
            for(int n = 0; n < colors.length; n++)
                if(color.equals(colors[n].getColor())) return colors[n];
        }
        else if(nTokens == 4)
        {
            float red = Float.parseFloat(stok.nextToken());
            float green = Float.parseFloat(stok.nextToken());
            float blue = Float.parseFloat(stok.nextToken());
            float alpha = Float.parseFloat(stok.nextToken());
            StsColor color = new StsColor(red, green, blue, alpha);
            for(int n = 0; n < colors.length; n++)
                if(color.equals(colors[n])) return colors[n];
        }
        else if(nTokens == 8)
        {
            stok.nextToken();
            float red = Float.parseFloat(stok.nextToken());
            stok.nextToken();
            float green = Float.parseFloat(stok.nextToken());
            stok.nextToken();
            float blue = Float.parseFloat(stok.nextToken());
            stok.nextToken();
            float alpha = Float.parseFloat(stok.nextToken());
            StsColor color = new StsColor(red, green, blue, alpha);
            for(int n = 0; n < colors.length; n++)
                if(color.equals(colors[n])) return colors[n];
        }
        StsException.systemError(StsColor.class, "colorFromString", "Failed parsing " + string);
        return StsColor.BLACK;
    }

    public boolean equals(StsColor color)
    {
        return StsMath.sameAs(red, color.red) && StsMath.sameAs(green, color.green) && StsMath.sameAs(blue, color.blue) && StsMath.sameAs(alpha, color.alpha);
    }
}





