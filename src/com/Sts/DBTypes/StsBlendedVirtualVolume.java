package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

public class StsBlendedVirtualVolume extends StsVirtualVolume implements StsTreeObjectI
{
    transient public static final int LT = 0;
    transient public static final int LE = 1;
    transient public static final int GT = 2;
    transient public static final int GE = 3;
    transient public static String[] LOGICALS =
        {"Less Than", "Less Than or Equal To", "Greater Than",
        "Greater Than or Equal To"};

    // these members are persistent, but not loaded from seis3d.txt.name file
    protected double condition = 1.0;
    protected int logical = UNDEFINED;
    static StsEditableColorscaleFieldBean vcolorscaleBean = new StsEditableColorscaleFieldBean(StsBlendedVirtualVolume.class,"colorscale");

    static public final StsFieldBean[] virtualDisplayFields =
        {
        new StsBooleanFieldBean(StsVirtualVolume.class, "isVisible", "Enable"),
        new StsBooleanFieldBean(StsVirtualVolume.class, "readoutEnabled", "Mouse Readout"),
        vcolorscaleBean
    };
    static public final StsFieldBean[] virtualPropertyFields = new StsFieldBean[]
        {
        new StsStringFieldBean(StsBlendedVirtualVolume.class, "name", true, "Name"),
        new StsStringFieldBean(StsBlendedVirtualVolume.class, "volumeOneName", false, "Seismic PostStack3d One"),
        new StsStringFieldBean(StsBlendedVirtualVolume.class, "volumeTwoName", false, "Seismic PostStack3d Two"),
        new StsDoubleFieldBean(StsBlendedVirtualVolume.class, "condition", false, "Condition"),
        new StsStringFieldBean(StsBlendedVirtualVolume.class, "logical", false, "Logical")
    };

    static StsObjectPanel virtualVolumeObjectPanel = null;
    
	public void setDataHistogram()
	{
		if (dataHist != null && colorscaleBean != null)
        {
            vcolorscaleBean.setHistogram(dataHist);
            vcolorscaleBean.revalidate();
        }
    }
	
    public StsBlendedVirtualVolume()
    {
        //System.out.println("BlendedVirtualVolume constructor called.");    	
    }
    
    public StsBlendedVirtualVolume(boolean persistent)
    {
        super(persistent);    	
    }
    
    public StsBlendedVirtualVolume(StsObject[] volumeList, String name, double condition, int logical)
    {
    	super(false);
    	
        StsToolkit.copyDifferentStsObjectAllFields(volumeList[0], this);
		clearNonRelevantMembers();

        setVolumes(volumeList);

        this.condition = condition;
        this.logical = logical;
        this.type = SEISMIC_BLEND;

        setName(name);
        colorscale = null;
        initializeColorscale();
        initialize(currentModel);
        isVisible = true;
        getBlendedVirtualVolumeClass().setIsVisibleOnCursor(true);
    }

    public float[] processBlock(float[] volOne, float[] volTwo, int nSamples)
    {
        try
        {
            float[] floatValues = new float[nSamples];
        	for (int n = 0; n < nSamples; n++)
        	{            
         	   if(volOne[n] == StsParameters.nullValue)
        	   {
        		   floatValues[n] = volTwo[n];
        		   continue;
        	   }  
         	   if(volTwo[n] == StsParameters.nullValue)
        	   {
        		   floatValues[n] = volOne[n];
        		   continue;
        	   }          	   
        		switch (logical)
        		{
        		case GT:
        			if (volTwo[n] > condition)
            			floatValues[n] = volOne[n];
            		else
            			floatValues[n] = volTwo[n];
        			break;
        		case GE:
            		if (volTwo[n] >= condition)
            			floatValues[n] = volOne[n];
            		else
            			floatValues[n] = volTwo[n];
            		break;
        		case LT:
        			if (volTwo[n] < condition)
        				floatValues[n] = volOne[n];
            		else
            			floatValues[n] = volTwo[n];
        			break;
        		case LE:
            		if (volTwo[n] <= condition)
            			floatValues[n] = volOne[n];
            		else
            			floatValues[n] = volTwo[n];
            		break;
        		default:
        			break;
        		}
        		if((floatValues[n] < dataMin) && (floatValues[n] != StsParameters.nullValue)) dataMin = floatValues[n];
        		if((floatValues[n] > dataMax) && (floatValues[n] != StsParameters.nullValue)) dataMax = floatValues[n];            
        	}
            return floatValues;        	
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.processBlock() failed.", e, StsException.WARNING);
            return null;
        }
    }
    
    public float[] processBlock(float[] volOne, int nSamples)
    {
        try
        {
           float[] floatValues = new float[nSamples];
           for (int n = 0; n < nSamples; n++)
           {
        	   if(volOne[n] == StsParameters.nullValue)
        	   {
        		   floatValues[n] = volOne[n];
        		   continue;
        	   } 	
                	
               switch (logical)
               {
                  case GT:
                	  if (volOne[n] > condition)
                		  floatValues[n] = volOne[n];
                	  else
                		  floatValues[n] = 0.0f;
                	  break;
                  case GE:
                	  if (volOne[n] >= condition)
                		  floatValues[n] = volOne[n];
                	  else
                		  floatValues[n] = 0.0f;
                	  break;
                  case LT:
                	  if (volOne[n] < condition)
                		  floatValues[n] = volOne[n];
                	  else
                		  floatValues[n] = 0.0f;
                	  break;
                  case LE:
                	  if (volOne[n] <= condition)
                		  floatValues[n] = volOne[n];
                	  else
                		  floatValues[n] = 0.0f;
                	  break;
                  default:
                	  floatValues[n] = 0.0f;
                  if((floatValues[n] < dataMin) && (floatValues[n] != StsParameters.nullValue)) dataMin = floatValues[n];
                  if((floatValues[n] > dataMax) && (floatValues[n] != StsParameters.nullValue)) dataMax = floatValues[n];
                }
           }
           return floatValues;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.processBlock() failed.", e, StsException.WARNING);
            return null;
        }
    }
    
    public byte[] processPlaneData(byte[][] planeValues)
    {
        try
        {
            int nPlanes = planeValues.length;
            int nPlanePoints = planeValues[0].length;

            StsVolumeDisplayable[] planeVolumes = new StsVolumeDisplayable[nPlanes];
            for (int n = 0; n < nPlanes; n++)
            {
                planeVolumes[n] = (StsVolumeDisplayable) volumes.getElement(n);
            }
            if (nPlanes == 1)
            {
                float[] floatValues = new float[nPlanePoints];
                int[] intValues = new int[nPlanePoints];

                StsVolumeDisplayable volume = planeVolumes[0];

                switch (logical)
                {
                    case GT:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                if (value > condition)
                                    floatValues[n] = value;
                                else
                                    floatValues[n] = 0.0f;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case GE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                if (value >= condition)
                                    floatValues[n] = value;
                                else
                                    floatValues[n] = 0.0f;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case LT:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                if (value < condition)
                                    floatValues[n] = value;
                                else
                                    floatValues[n] = 0.0f;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case LE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                if (value <= condition)
                                    floatValues[n] = value;
                                else
                                    floatValues[n] = 0.0f;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    default:
                        return planeValues[0];
                }
                return computeUnscaledByteValues(floatValues);

            }
            else if(nPlanes == 2)
            {
                StsVolumeDisplayable volume0 = planeVolumes[0];
                StsVolumeDisplayable volume1 = planeVolumes[1];
                float[] floatValues = new float[nPlanePoints];

                switch (logical)
                {
                    case GT:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            float value0 = volume0.getScaledValue(planeValues[0][n]);
                            float value1 = volume1.getScaledValue(planeValues[1][n]);
                            if (value1 > condition)
                                floatValues[n] = value0;
                            else
                                floatValues[n] = value1;
                        }
                        break;
                    case GE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            float value0 = volume0.getScaledValue(planeValues[0][n]);
                            float value1 = volume1.getScaledValue(planeValues[1][n]);
                            if (value1 >= condition)
                                floatValues[n] = value0;
                            else
                                floatValues[n] = value1;
                        }
                        break;
                    case LT:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            float value0 = volume0.getScaledValue(planeValues[0][n]);
                            float value1 = volume1.getScaledValue(planeValues[1][n]);
                            if (value1 < condition)
                                floatValues[n] = value0;
                            else
                                floatValues[n] = value1;
                        }
                        break;
                    case LE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            float value0 = volume0.getScaledValue(planeValues[0][n]);
                            float value1 = volume1.getScaledValue(planeValues[1][n]);
                            if (value1 <= condition)
                                floatValues[n] = value0;
                            else
                                floatValues[n] = value1;
                        }
                        break;
                    default:
                        return planeValues[0];
                }
                return computeUnscaledByteValues(floatValues);
            }
            else
            {
                return planeValues[0];
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.processBytePlaneData() failed.",
                                         e, StsException.WARNING);
            return null;
        }
    }

    public String getLogical()
    {
        String operString = null;
        switch (logical)
        {
            case LT:
                operString = "Less Than";
                break;
            case LE:
                operString = "Less Than or Equal To";
                break;
            case GT:
                operString = "Greater Than";
                break;
            case GE:
                operString = "Greater Than or Equal To";
                break;
            default:
                operString = "Undefined";
        }
        return operString;
    }

    public void setLogical(String operString)
    {
        if (operString.equals("Less Than") || operString.equals("LT"))
        {
            logical = LT;
        }
        else if (operString.equals("Less Than or Equal To") || operString.equals("LE"))
        {
            logical = LE;
        }
        else if (operString.equals("Greater Than") || operString.equals("GT"))
        {
            logical = GT;
        }
        else if (operString.equals("Greater Than or Equal To") || operString.equals("GE"))
        {
            logical = GE;
        }
        else
        {
            logical = UNDEFINED;
        }
        return;
    }

    public double getCondition()
    {return condition;
    }

    public void setCondition(double cond)
    {
        condition = cond;
        return;
    }
    public StsFieldBean[] getDisplayFields()
    {
        return virtualDisplayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return virtualPropertyFields;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (virtualVolumeObjectPanel == null)
        {
            virtualVolumeObjectPanel = StsObjectPanel.constructor(this, true);
        }
        return virtualVolumeObjectPanel;
    }

    static public StsBlendedVirtualVolumeClass getBlendedVirtualVolumeClass()
    {
        return (StsBlendedVirtualVolumeClass) currentModel.getCreateStsClass(StsBlendedVirtualVolume.class);
    }
}