//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.UI.Beans.*;
import com.Sts.UI.Sounds.*;
import com.Sts.Utilities.*;

public class StsPump extends StsStaticSensor
{
    protected byte pumpType = TYPE1;
    protected float maxPumpVolume = 1000.0f;
    protected float warningRate = 75.0f;
    protected float criticalRate = 90.0f;
    protected boolean pumpOn = true;
    
    // Event symbols
    public static final byte TYPE1 = 0;
    public static final byte TYPE2 = 1;
    public static final byte TYPE3 = 2;
    public static final byte TYPE4 = 3;
    static public final String[] PUMP_TYPE_STRINGS = new String[] { "Type1", "Type2", "Type3", "Type4"};
    static public final byte[] PUMP_TYPES = new byte[] { TYPE1, TYPE2, TYPE3, TYPE4 };
    
    static protected StsComboBoxFieldBean pumpTypeListBean;
    static protected StsDoubleFieldBean pumpXBean;
    static protected StsDoubleFieldBean pumpYBean;
    static protected StsDoubleFieldBean pumpZBean;
    static protected StsComboBoxFieldBean pumpDisplayTypeBean;    
    /** default constructor */
    public StsPump()
    {
    	super();
    }

    public StsPump(StsWell well, String name)
    {
        this(well, name, 0.0, 0.0, 0.0);
    }

    public StsPump(StsWell well, String name, double x, double y, double z)
    {
    	super(well, name, x, y, z);
    	setSymbolString(SYMBOL_TYPE_STRINGS[SPHERE]);
        scaleMin = 0.0f;
    	scaleMax = 100.0f;
    }

    static public StsPump nullSensorConstructor(String name)
    {
        return new StsPump(null, name);
    }
    
    public StsFieldBean[] getDisplayFields()
    {
       try
       {
           if (displayFields == null)
           {
               pumpTypeListBean = new StsComboBoxFieldBean(StsPump.class, "pumpTypeString", "Type:", PUMP_TYPE_STRINGS);
               pumpXBean = new StsDoubleFieldBean(StsPump.class, "xLoc", true, "X:", false);
               pumpYBean = new StsDoubleFieldBean(StsPump.class, "yLoc", true, "Y:", false);
               pumpZBean = new StsDoubleFieldBean(StsPump.class, "zLoc", true, "Depth:", false);
               displayFields = new StsFieldBean[]
               {
                    new StsBooleanFieldBean(StsPump.class, "isVisible", "Enable"),
                    pumpTypeListBean,
            		new StsFloatFieldBean(StsPump.class, "warningRate", true, "Warn Volume"),
         		    new StsFloatFieldBean(StsPump.class, "criticalRate", true, "Critical Volume"),
            		new StsFloatFieldBean(StsPump.class, "maxPumpVolume", true, "Maximum Volume (cuft/sec)"),
         		    pumpDisplayTypeBean,
                    pumpXBean, pumpYBean, pumpZBean
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
    public float getMaxPumpVolume() { return maxPumpVolume; }
    public float getWarningRate() { return warningRate; }
    public float getCriticalRate() { return criticalRate; }
    
    public void setMaxPumpVolume(float volume) { maxPumpVolume = volume;         dbFieldChanged("maxPumpVolume", maxPumpVolume);}
    public void setWarningRate(float rate) { warningRate = rate;         dbFieldChanged("warningRate", warningRate);}
    public void setCriticalRate(float rate) { criticalRate = rate;         dbFieldChanged("criticalRate", criticalRate);}
    
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
    	StsPumpClass pumpClass = (StsPumpClass)currentModel.getCreateStsClass(StsPump.class);
    	if((value > warningRate) && (value < criticalRate))
    	{
    		color = StsColor.YELLOW;
    		if(pumpClass.getEnableSound())
    			StsSound.play(StsSound.BEEP1);
    	}
    	else if(value > criticalRate)
    	{
    		color = StsColor.RED;
    		if(pumpClass.getEnableSound())
    			StsSound.play(StsSound.BEEP4);
    	}	
    	return color;
    }   

    public String getPumpTypeString()
    {
        return PUMP_TYPE_STRINGS[pumpType];
    }

    public void setPumpTypeString(String stype)
    {
        for(int i=0; i<PUMP_TYPE_STRINGS.length; i++)
        {
            if(stype == PUMP_TYPE_STRINGS[i])
            {
                if(pumpType == i) return;
                pumpType = (byte)i;
                dbFieldChanged("pumpType", pumpType);
                currentModel.viewObjectRepaint(this, this);
                return;
            }
        }
    } 
    public boolean getEnablePump()
    {
    	return pumpOn;
    }    
    public void enablePump(boolean val)
    {
    	pumpOn = val;
    }

}