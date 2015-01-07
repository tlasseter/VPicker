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

public class StsMathVirtualVolume extends StsVirtualVolume implements StsTreeObjectI
{
    transient public static final int UNDEFINED = 0;
    transient public static final int ADDITION = 1;
    transient public static final int SUBTRACTION = 2;
    transient public static final int MULTIPLY = 3;
    transient public static final int DIVIDE = 4;
    transient public static final int AVERAGE = 5;
    transient public static final int MAXIMUM = 6;
    transient public static final int MINIMUM = 7;
    transient public static final int ABSOLUTE = 8;
    transient public static final int POWER = 9;
    transient public static final int LOG = 10;
    transient public static final int ASQRT = 11;
    transient public static final int SQUARE = 12;
    transient public static final int MODULUS = 13;
    transient public static final int INVERSE = 14;
    transient public static final int SIN = 15;
    transient public static final int COS = 16;
    transient public static final int TRANSPARENCY = 17;
    transient public static final int SCALAR = 18;
    transient public static String[] OPERATORS =
        {"undefined", "+ (add)", "- (subtract)", "* (multiply)", "/ (divide)",
        "average", "maximum", "minimum", "absolute", "power of", "log", "absolute sqrt",
        "square", "modulus", "inverse", "sin", "cos", "transparency", "scalar"};

    // these members are persistent, but not loaded from seis3d.txt.name file
    protected double scalar = 1.0;
    protected int operator = UNDEFINED;
    protected double condition = 1.0;

    static StsObjectPanel virtualVolumeObjectPanel = null;
    static StsEditableColorscaleFieldBean vcolorscaleBean = new StsEditableColorscaleFieldBean(StsMathVirtualVolume.class,"colorscale");

    static public final StsFieldBean[] virtualDisplayFields =
        {
        new StsBooleanFieldBean(StsVirtualVolume.class, "isVisible", "Enable"),
        new StsBooleanFieldBean(StsVirtualVolume.class, "readoutEnabled", "Mouse Readout"),
                vcolorscaleBean
    };
    static public final StsFieldBean[] virtualPropertyFields = new StsFieldBean[]
        {
        new StsStringFieldBean(StsMathVirtualVolume.class, "name", true, "Name"),
        new StsStringFieldBean(StsMathVirtualVolume.class, "volumeOneName", false, "Seismic PostStack3d One"),
        new StsStringFieldBean(StsMathVirtualVolume.class, "volumeTwoName", false, "Seismic PostStack3d Two"),
        new StsDoubleFieldBean(StsMathVirtualVolume.class, "scalar", false, "Scalar"),
        new StsStringFieldBean(StsMathVirtualVolume.class, "operator", false, "Operator"),
    };

    public StsMathVirtualVolume()
    {
        //System.out.println("MathVirtualVolume constructor called.");
    }
    
	public void setDataHistogram()
	{
		if (dataHist != null && colorscaleBean != null)
        {
            vcolorscaleBean.setHistogram(dataHist);
            vcolorscaleBean.revalidate();
        }
    }
	
	public StsMathVirtualVolume(boolean persistent)
	{
		super(persistent);
	}

    public StsMathVirtualVolume(StsObject[] volumeList, String name, int oper, double scalar)
    {
        super(false);

        StsToolkit.copyDifferentStsObjectAllFields(volumeList[0], this);
        clearNonRelevantMembers();

        setVolumes(volumeList);

        this.scalar = scalar;
        this.operator = oper;
        this.type = SEISMIC_MATH;

        setName(name);
        colorscale = null;
        initializeColorscale();
        initialize(currentModel);
        //initializeData();
        isVisible = true;
        getMathVirtualVolumeClass().setIsVisibleOnCursor(true);
    }

    public float[] processBlock(float[] valOne, int size)
    {
    	float[] results = new float[size];
    	for(int i=0; i<size; i++)
    	{
            if(valOne[i] == StsParameters.nullValue) 
            	results[i] = valOne[i];
            
            switch (operator)
            {
                case SCALAR:
                    results[i] = (float)(scalar * valOne[i]);
                    break;
                case ABSOLUTE:
                	results[i] =  Math.abs(valOne[i]);
                	break;
                case LOG:
                	results[i] =  (float)Math.log(valOne[i]);
                	break;
                case ASQRT:
                	results[i] =  (float)Math.sqrt(Math.abs(valOne[i]));
                	break;
                case SQUARE:
                	results[i] =  valOne[i] * valOne[i];
                	break;
                case INVERSE:
                    if (valOne[i] == 0.0f)
                    	results[i] =  0.0f;
                    else
                    	results[i] =  1.0f / valOne[i];
                	break;
                case SIN:
                	results[i] =  (float)Math.sin(valOne[i]);
                	break;
                case COS:
                	results[i] =  (float)Math.cos(valOne[i]);
                	break;
                case ADDITION:
                	results[i] =  (float)(valOne[i] + scalar);
                	break;
                case SUBTRACTION:
                	results[i] =  (float)(valOne[i] - scalar);
                	break;
                case MULTIPLY:
                	results[i] =  (float)(valOne[i] * scalar);
                	break;
                case DIVIDE:
                    if (scalar == 0.0)
                    	results[i] =  StsParameters.nullValue;
                    else
                    	results[i] =  (float)(valOne[i] / scalar);
                	break;
                case AVERAGE:
                	results[i] =  (float)((valOne[i] - scalar)) / 2.0f;
                	break;
                case MAXIMUM:
                	results[i] =  (float)Math.max(valOne[i], scalar);
                	break;
                case MINIMUM:
                	results[i] =  (float)Math.min(valOne[i], scalar);
                	break;
                case POWER:
                	results[i] =  (float)Math.pow(valOne[i], scalar);
                	break;
                case MODULUS:
                	results[i] =  (float)(valOne[i] % scalar);
                	break;
                case TRANSPARENCY:
                	results[i] =  0.0f;
                	break;
                default:
                	results[i] =  0.0f;
                	break;
            }
            if((results[i] < dataMin) && (results[i] != StsParameters.nullValue)) dataMin = results[i];
            if((results[i] > dataMax) && (results[i] != StsParameters.nullValue)) dataMax = results[i];            
    	}
    	return results;
    }    

    public float[] processBlock(float[] valOne, float[] valTwo, int size)
    {
    	float[] results = new float[size];
    	for(int i=0; i<size; i++)
    	{
    		if((valOne[i] == StsParameters.nullValue) || (valTwo[i] == StsParameters.nullValue))
    		{
    			results[i] = StsParameters.nullValue;
    			continue;
    		}
    		switch (operator)
    		{
    		case ADDITION:
    			results[i] = valOne[i] + valTwo[i];
           	   	break;
    		case SUBTRACTION:
    			results[i] = valOne[i] - valTwo[i];
    			break;
    		case MULTIPLY:
    			results[i] = valOne[i] * valTwo[i];
           	   	break;
    		case DIVIDE:
    			if (valTwo[i] == 0)
    				results[i] = StsParameters.nullValue;
    			else
    				results[i] = valOne[i] / valTwo[i];
    			break;
    		case AVERAGE:
    			results[i] = (valOne[i] - valTwo[i]) / 2.0f;
           	   	break;
    		case MAXIMUM:
    			results[i] = Math.max(valOne[i], valTwo[i]);
           	   	break;
    		case MINIMUM:
    			results[i] = Math.min(valOne[i], valTwo[i]);
           	   	break;
    		case POWER:
    			results[i] = (float) Math.pow(valOne[i], valTwo[i]);
           	   	break;
    		case MODULUS:
    			results[i] = valOne[i] % valTwo[i];
    			break;
    		case UNDEFINED:
    			results[i] = 0.0f;
           	   	break;
    		default:
    			results[i] = 0.0f;
       	   	   	break;
    		}  
    		if((results[i] < dataMin) && (results[i] != StsParameters.nullValue)) dataMin = results[i];
    		if((results[i] > dataMax) && (results[i] != StsParameters.nullValue)) dataMax = results[i];
    	}
    	return results;
    }
    
    public byte[] processPlaneData(byte[][] planeValues)
    {
        try
        {
            if(planeValues == null)
                return null;
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

                switch (operator)
                {
                    case SCALAR:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                                floatValues[n] = (float) scalar * volume.getScaledValue(planeValues[0][n]);
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case ABSOLUTE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = Math.abs(value);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case LOG:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                if(value <= 0.0)
                                    floatValues[n] = StsParameters.nullValue;
                                else
                                    floatValues[n] = (float) Math.log(value);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case ASQRT:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = (float) Math.sqrt(Math.abs(value));
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case SQUARE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = value * value;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case INVERSE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                if (value == 0.0f)
                                    floatValues[n] = 0.0f;
                                else
                                    floatValues[n] = 1.0f / value;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case SIN:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = (float)Math.sin((double)value);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case COS:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = (float)Math.cos((double)value);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case ADDITION:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = value0 + (float)scalar;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case SUBTRACTION:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = value0 - (float)scalar;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MULTIPLY:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = value0 * (float)scalar;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case DIVIDE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                if (scalar == 0.0)
                                    floatValues[n] = StsParameters.nullValue;
                                else
                                    floatValues[n] = value0 / (float)scalar;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case AVERAGE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = (value0 - (float)scalar) / 2.0f;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MAXIMUM:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = Math.max(value0, (float)scalar);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MINIMUM:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = Math.min(value0, (float)scalar);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case POWER:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = (float) Math.pow(value0, (float)scalar);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MODULUS:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(planeValues[0][n] != -1)
                            {
                                float value0 = volume.getScaledValue(planeValues[0][n]);
                                floatValues[n] = value0 % (float)scalar;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case TRANSPARENCY:
                        return planeValues[0];
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

                switch (operator)
                {
                    case ADDITION:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = value0 + value1;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case SUBTRACTION:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = value0 - value1;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MULTIPLY:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = value0 * value1;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case DIVIDE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                if (value1 == 0)
                                    floatValues[n] = StsParameters.nullValue;
                                else
                                    floatValues[n] = value0 / value1;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case AVERAGE:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = (value0 - value1) / 2.0f;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MAXIMUM:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = Math.max(value0, value1);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MINIMUM:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = Math.min(value0, value1);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case POWER:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = (float) Math.pow(value0, value1);
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case MODULUS:
                        for (int n = 0; n < nPlanePoints; n++)
                        {
                            if(!isAValueNull(planeValues[0][n], planeValues[1][n]))
                            {
                                float value0 = volume0.getScaledValue(planeValues[0][n]);
                                float value1 = volume1.getScaledValue(planeValues[1][n]);
                                floatValues[n] = value0 % value1;
                            }
                            else
                                floatValues[n] = StsParameters.nullValue;
                        }
                        break;
                    case UNDEFINED:
                    default:
                        return planeValues[0];
                }
                return computeUnscaledByteValues(floatValues);
            }
            else
                return planeValues[0];
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.processBytePlaneData() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public String getOperator()
    {
        String operString = null;
        switch (operator)
        {
            case ADDITION:
                operString = "Add";
                break;
            case SUBTRACTION:
                operString = "Subtract";
                break;
            case MULTIPLY:
                operString = "Multiply";
                break;
            case DIVIDE:
                operString = "Divide";
                break;
            case AVERAGE:
                operString = "Average";
                break;
            case MAXIMUM:
                operString = "Maximum";
                break;
            case MINIMUM:
                operString = "Minimum";
                break;
            case ABSOLUTE:
                operString = "Absolute";
                break;
            case POWER:
                operString = "Power Of";
                break;
            case LOG:
                operString = "Log";
                break;
            case ASQRT:
                operString = "Absolute Sqrt";
                break;
            case SQUARE:
                operString = "Square";
                break;
            case MODULUS:
                operString = "Modulus";
                break;
            case INVERSE:
                operString = "Inverse";
                break;
            case SIN:
                operString = "Sin";
                break;
            case COS:
                operString = "Cos";
                break;
            case TRANSPARENCY:
                operString = "Transparency";
                break;
            case SCALAR:
                operString = "Scalar";
                break;
            case UNDEFINED:
            default:
                operString = "Undefined";
        }
        return operString;
    }

    public void setOperator(String operString)
    {
        if (operString.equals("Add"))
        {
            operator = ADDITION;
        }
        else if (operString.equals("Subtract"))
        {
            operator = SUBTRACTION;
        }
        else if (operString.equals("Multiply"))
        {
            operator = MULTIPLY;
        }
        else if (operString.equals("Divide"))
        {
            operator = DIVIDE;
        }
        else if (operString.equals("Average"))
        {
            operator = AVERAGE;
        }
        else if (operString.equals("Maximum"))
        {
            operator = MAXIMUM;
        }
        else if (operString.equals("Minimum"))
        {
            operator = MINIMUM;
        }
        else if (operString.equals("Absolute"))
        {
            operator = ABSOLUTE;
        }
        else if (operString.equals("Power Of"))
        {
            operator = POWER;
        }
        else if (operString.equals("Log"))
        {
            operator = LOG;
        }
        else if (operString.equals("Absolute Sqrt"))
        {
            operator = ASQRT;
        }
        else if (operString.equals("Square"))
        {
            operator = SQUARE;
        }
        else if (operString.equals("Modulus"))
        {
            operator = MODULUS;
        }
        else if (operString.equals("Inverse"))
        {
            operator = INVERSE;
        }
        else if (operString.equals("Sin"))
        {
            operator = SIN;
        }
        else if (operString.equals("Cos"))
        {
            operator = COS;
        }
        else if (operString.equals("Transparency"))
        {
            operator = TRANSPARENCY;
        }
        else if (operString.equals("Scalar"))
        {
            operator = SCALAR;
        }
        else
        {
            operator = UNDEFINED;
        }
        return;
    }

    public double getScalar()
    {return scalar;
    }

    public void setScalar(double scale)
    {
        scalar = scale;
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

    static public StsMathVirtualVolumeClass getMathVirtualVolumeClass()
    {
        return (StsMathVirtualVolumeClass) currentModel.getCreateStsClass(StsMathVirtualVolume.class);
    }
}