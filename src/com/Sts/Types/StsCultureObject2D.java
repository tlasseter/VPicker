package com.Sts.Types;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;

public abstract class StsCultureObject2D extends StsMainTimeObject
{
    public StsColor stsColor = new StsColor(StsColor.RED);
    public int size = 2;
    public int decimation = 10;
    public int scale = 2;
    public int width = 2;
    public byte symbolType = PT2D;
    public StsColorscale colorscale;
    public boolean planar = true;

    public static final byte D2 = 0;
    public static final byte D3 = 1;

    public static final byte PT2D = 0;
    public static final byte SPHERE = 1;
    public static final byte CUBE = 2;
    public static final byte CYLINDER = 3;
    public static String[] typeStrings = { "2D", "SPHERE", "CUBE", "CYLINDER" };
    static public final byte[] SYMBOL_TYPES = new byte[] { PT2D, SPHERE, CUBE, CYLINDER };

    public static final byte LINE = 0;
    public static final byte POINT = 1;
    public static final byte XY = 2;
    public static final byte TEXT = 3;
    transient public static String[] TYPES = {"LINE", "POINT", "XY", "TEXT"};

    abstract public void draw(StsCultureDisplayable cultureDisplayable, boolean mapToSurface, byte origZDomain, StsGLPanel3d glPanel3d);
    abstract public void draw2d(StsCultureDisplayable cultureDisplayable, boolean mapToSurface, byte origZDomain, StsGLPanel glPanel);
    public void initialize() {}
    public boolean initialize(StsModel model){ return true; }
    abstract public float getMaxDepth(byte zDomainOriginal);
    abstract public float getMaxTime(byte zDomainOriginal);
    abstract public float getMinDepth(byte zDomainOriginal);
    abstract public float getMinTime(byte zDomainOriginal);

    public void setPlanar(boolean value) 
    { 
    	planar = value;
		fieldChanged("planar", planar);    	
    }
    public boolean isPlanar() { return planar; }
    
    abstract public double getXMin();
    abstract public double getYMin();
    abstract public double getXMax();
    abstract public double getYMax();
    abstract public StsPoint getPointAt(int index);
    abstract public int getNumPoints();

    public StsColor getStsColor() { return stsColor; }
    public void setStsColor(StsColor color)
    {
        this.stsColor = color;
		fieldChanged("stsColor", color);        
    }
    public int getSize() { return size; }
    public void setSize(int size) 
    { 
    	fieldChanged("size", size); 
    	this.size = size; 
    }
    public void setScale(int scale) 
    { 
    	fieldChanged("scale", scale);
    	this.scale = scale; 
    }
    public int getScale() { return scale; }
    public void setWidth(int scale) 
    { 
    	fieldChanged("width", scale); 
    	this.width = scale; 
    }
    public int getWidth() { return width; }
    public void setSymbolType(byte type) { this.symbolType = type; }
    public byte getSymbolType() { return symbolType; }
    public void setSymbolTypeByString(String typeString)
    {
        for(int i=0; i<typeStrings.length; i++)
        {
            if(typeStrings[i].equalsIgnoreCase(typeString))
            {
                symbolType = (byte)i;
        		fieldChanged("symbolType", symbolType);                
                return;
            }
        }
    }
    public String[] getTypeStrings() { return typeStrings; }
}
