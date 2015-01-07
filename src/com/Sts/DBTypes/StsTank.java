//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;


import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Sounds.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsTank extends StsStaticSensor
{
    protected byte tankType = TYPE1;
    protected float tankRadius = 10.0f;
    protected float tankHeight = 100.0f;
    protected float warningVolume = 75.0f;
    protected float criticalVolume = 95.0f;
    protected boolean onZSlice = false;
    protected StsPump pump = null;

    // Event symbols
    public static final byte TYPE1 = 0;
    public static final byte TYPE2 = 1;
    public static final byte TYPE3 = 2;
    public static final byte TYPE4 = 3;
    static public final String[] TANK_TYPE_STRINGS = new String[] { "Cylindrical", "Spherical", "Cubic", "Unknown"};
    static public final byte[] TANK_TYPES = new byte[] { TYPE1, TYPE2, TYPE3, TYPE4 };
    
    static protected StsComboBoxFieldBean tankTypeListBean;
    static protected StsDoubleFieldBean tankXBean;
    static protected StsDoubleFieldBean tankYBean;
    static protected StsDoubleFieldBean tankZBean;
    static protected StsComboBoxFieldBean tankDisplayTypeBean;
    
    /** default constructor */
    public StsTank()
    {
    	super();
    }

    public StsTank(StsWell well, String name, boolean persistent)
    {
        super(well, name);
    	setDisplayTypeString(displayTypeStrings[SHOW_SINGLE]);
    	setSymbolString(SYMBOL_TYPE_STRINGS[CYLINDER]);
    	scaleMin = 0.0f;
    	scaleMax = 100.0f;        
    }

    public StsTank(StsWell well, String name)
    {
        this(well, name, 0.0, 0.0, 0.0);
    }

    public StsTank(StsWell well, String name, double x, double y, double z)
    {
    	super(well, name, x, y, z);
    }

    static public StsTank nullSensorConstructor(String name)
    {
        return new StsTank(null, name, false);
    }
    
    public StsFieldBean[] getDisplayFields()
    {
       try
       {
           if (displayFields == null)
           {
               tankTypeListBean = new StsComboBoxFieldBean(StsTank.class, "tankTypeString", "Type:", TANK_TYPE_STRINGS);
               tankXBean = new StsDoubleFieldBean(StsTank.class, "xLoc", true, "X:", false);
               tankYBean = new StsDoubleFieldBean(StsTank.class, "yLoc", true, "Y:", false);
               tankZBean = new StsDoubleFieldBean(StsTank.class, "zLoc", true, "Depth:", false);
               tankDisplayTypeBean = new StsComboBoxFieldBean(StsTank.class, "displayTypeString", "Show:", displayTypeStrings);               
               displayFields = new StsFieldBean[]
               {
                    new StsBooleanFieldBean(StsTank.class, "isVisible", "Enable"),
                    new StsBooleanFieldBean(StsTank.class, "onZSlice", "Plot on Slice"),
                    tankTypeListBean,
            		new StsFloatFieldBean(StsTank.class, "warningVolume", true, "Warn Volume"),
         		    new StsFloatFieldBean(StsTank.class, "criticalVolume", true, "Critical Volume"),
            		new StsFloatFieldBean(StsTank.class, "tankHeight", true, "Tank Height"),
         		    new StsFloatFieldBean(StsTank.class, "tankRadius", true, "Tank Radius"),
         		    tankDisplayTypeBean,
                    tankXBean, tankYBean, tankZBean
               };
           }
            return displayFields;
        }
        catch (Exception e)
        {
            StsException.outputException("StsTank.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    /** add to reference lists */
    public void addTimeCurves(StsTimeCurve[] timeCurves)
    {
        if(timeCurves == null)
            return;

        for (int n = 0; n < timeCurves.length; n++)
            addTimeCurve(timeCurves[n]);

    	initialize();

        dataMin = timeCurves[firstNonPositionalCurve()].getCurveMin();
        dataMax = timeCurves[firstNonPositionalCurve()].getCurveMax();

        if(getTimeCurve(StsLogVector.X) != null)
        {
            xOrigin = getTimeCurve(StsLogVector.X).getValueVector().getOrigin();
            yOrigin = getTimeCurve(StsLogVector.Y).getValueVector().getOrigin();
        }
        else
        {
            xOrigin = 0.0f;
            yOrigin = 0.0f;
        }
    }   
    public boolean getOnZSlice() { return onZSlice; }
    public float getTankHeight() { return tankHeight; }
    public float getTankRadius() { return tankRadius; }
    public float getWarningVolume() { return warningVolume; }
    public float getCriticalVolume() { return criticalVolume; }
    
    public void setOnZSlice(boolean val) { onZSlice = val;         dbFieldChanged("onZSlice", onZSlice);} 
    public void setTankHeight(float height) { tankHeight = height;         dbFieldChanged("tankHeight", tankHeight);}
    public void setTankRadius(float radius) { tankRadius = radius;         dbFieldChanged("tankRadius", tankRadius);}
    public void setWarningVolume(float warning) { warningVolume = warning;         dbFieldChanged("warningVolume", warningVolume);}
    public void setCriticalVolume(float critical) { criticalVolume = critical;         dbFieldChanged("criticalVolume", criticalVolume);}
    public void setCoordinateType(byte type)
    {
        // Verify that the X, Y and Z vectors exist
        if(type == DYNAMIC)
        {
            if(!canBeDynamic())
                type = NONE;
        }
        coordinateType = type;
        if((type == DYNAMIC) || (type == NONE))
        {
            tankXBean.setEditable(false);
            tankYBean.setEditable(false);
            tankZBean.setEditable(false);
        }
        else
        {
            tankXBean.setEditable(true);
            tankYBean.setEditable(true);
            tankZBean.setEditable(true);
        }
    }    
    public String getTankTypeString()
    {
        return TANK_TYPE_STRINGS[tankType];
    }

    public void setTankTypeString(String stype)
    {
        for(int i=0; i<TANK_TYPE_STRINGS.length; i++)
        {
            if(stype == TANK_TYPE_STRINGS[i])
            {
                if(tankType == i) return;
                tankType = (byte)i;
                dbFieldChanged("tankType", tankType);
                currentModel.viewObjectRepaint(this, this);
                return;
            }
        }
    }
    
    public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean pick)
    {
        if(!verifyDisplayable()) return;
        
        GL gl = glPanel3d.getGL();

        // Temporarily ignore settings for duration and type if user is interactively selecting a region on a time series plot.
        long duration = this.displayDuration;
        byte type = this.displayType;

        if(!checkLoadVectors())
        	return;

        int points = 0;
        try
        {
            float[] xyz;
            float[] xy;
            glu = glPanel3d.getGLU();
            float sizeFloat, colorFloat;
            
            for (int i = 0; i < timeLongs.length; i++)
            {
                if(sizeFloats == null)
                    sizeFloat = getNumberBins();
                else
                    sizeFloat = sizeFloats[i];
                
                if(colorFloats == null)
                    colorFloat = 1;
                else
                    colorFloat = colorFloats[i];

                // Determine the xyz values.
               float zValue = (float)getZLoc() + zShift;
                if(onZSlice)
                	zValue = currentModel.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
               float xValue = (float)getXLoc() + xShift;
               float yValue = (float)getYLoc() + yShift;
               xyz = getXYZValue(xValue, yValue, zValue);

               // If time is enabled, verify the the time of the event is before the current time and in range.
               if(getStaticSensorClass().getEnableTime())
               {
                   // Check if before current time
                   if(!checkTime(timeLongs, i, points, duration, type) && !clustering)
               	        continue;
                   // Check if in range
                   if(!checkValue(timeLongs[i]) && !clustering)
               	        continue;
               }
               points++;

                // Draw the point
                if(sizeFloats != null)
                    sizeFloat = computeSize(scaleMin, scaleMax, sizeFloats[i]);
                if(pick) gl.glPushName(i);

                final StsColor color = defineColor(sizeFloat);
                drawTank(glPanel3d, xyz, sizeFloat, color);

                // Draw Name
                xyz[2] = xyz[2] - sizeFloat;                
                if(displayName) StsGLDraw.fontHelvetica12(gl, xyz, name);
                
                if(pick) gl.glPopName();
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSensor.display failed.", e, StsException.WARNING);
            return;
        }
    } 

    /** draw a solid color from the base up to the level height; and transparent from there to the top of the tank.
     *  Positive is down, so we flip height values to negative.
     */
    public void drawTank(StsGLPanel3d glPanel3d, float[] baseXyz, float levelHeight, StsColor color)
    {
    	float zscale = glPanel3d.getZScale();
        StsGLDraw.drawCylinder(glPanel3d, baseXyz, color, 50, zscale*levelHeight);
        float[] levelXyz = new float[] { baseXyz[0], baseXyz[1], baseXyz[2]-levelHeight };
        StsColor transparentColor = new StsColor(0.5f, 0.5f, 0.5f, 0.5f);
        float emptyHeight = 100.0f - levelHeight;
        StsGLDraw.drawCylinder(glPanel3d, levelXyz, transparentColor, 50, zscale*emptyHeight);
    }
    
    public void setProperty(StsTimeCurve vector)
    {
        vector = findPropertyVectorByName(vector);
        if((vector == null) || (propertyVector == vector))
            return;
        propertyVector = vector;
        propertyVectorName = vector.getName();
        dbFieldChanged("propertyVectorName", propertyVectorName);
        currentModel.viewObjectRepaint(this, this);
    }
    
    public StsColor defineColor(float value)
    {
    	StsColor color = StsColor.GREEN;
    	StsTankClass tankClass = (StsTankClass)currentModel.getCreateStsClass(StsTank.class);
    	if((value > warningVolume) && (value < criticalVolume))
    	{
    		color = StsColor.YELLOW;
    		if(pump != null) 
    		{
    			if(pump.getEnablePump())
    	    		StsMessageFiles.infoMessage("Tank (" + this.getName() + ") is nearing full, reduce pump rate for " + pump.getName());
    		}    		
    		if(tankClass.getEnableSound())
    			StsSound.play(StsSound.BEEP1);
    	}
    	else if(value > criticalVolume)
    	{
    		if(pump != null) 
    		{
    			if(pump.getEnablePump())
    			{
    	    		StsMessageFiles.errorMessage("Tank (" + this.getName() + ") is full, shutting down pump " + pump.getName());    		    				
    				pump.enablePump(false);
    			}
    		}
    		color = StsColor.RED;
    		if(tankClass.getEnableSound())
    			StsSound.play(StsSound.BEEP4);
    	}
    	else
    	{
    		if(pump != null) 
    		{
    			if(!pump.getEnablePump())
    			{
    				StsMessageFiles.infoMessage("Tank (" + this.getName() + ") has capacity, activating pump " + pump.getName());    		   			
    				pump.enablePump(true);
    			}
    		}
    	}
    	return color;
    } 
    
    public StsPump getPump() { return pump; }
    public void setPump(StsPump pump)
    {
    	this.pump = pump;
        dbFieldChanged("pump", pump);    	
    }
}