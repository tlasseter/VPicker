
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.DB.DBCommand.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.lang.reflect.*;


public class StsSerialize implements StsSerializable, Serializable
{
    public StsSerialize()
    {
    }

    /** Used by all StsObject subTypes to access currentModel */
	protected static StsModel currentModel;
	public static boolean isDepth = true;

   	static public void setCurrentModel(StsModel model) { currentModel = model; }
   	static public StsModel getCurrentModel() { return currentModel; }

    static public Frame getFrame() 
    {
        if(currentModel == null) return null;
        return currentModel.win3d;
    }

    static public void setIsDepth(boolean b) { isDepth = b; }
	static public boolean getIsDepth() { return isDepth; }

    static public StsProject getCurrentProject()
    {
        if (currentModel==null) return null;
        return currentModel.getProject();
    }

    static public byte getCurrentZDomain()
    {
        if (currentModel==null) return StsProject.TD_NONE;
        return currentModel.getProject().getZDomain();
    }

	public void outputWarningException(String methodName, String message, Exception e)
	{
		outputException(methodName, message, e, StsException.WARNING);
	}

	public void outputFatalException(String methodName, String message, Exception e)
	{
		outputException(methodName, message, e, StsException.FATAL);
	}

	public void outputException(String methodName, String message, Exception e, int level)
	{
		StsException.outputException(getClass().getName() + "." + methodName + " " + message, e, level);
	}

	public void dbFieldChanged(String fieldName, boolean value)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Boolean(value), false);
	}

	public void dbFieldChanged(String fieldName, int value)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Integer(value), false);
	}

	public void dbFieldChanged(String fieldName, byte value)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Byte(value), false);
	}

	public void dbFieldChanged(String fieldName, char value)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Character(value), false);
	}

    public void dbFieldChanged(String fieldName, double value)
    {
        if (currentModel == null)
        {
            return;
        }
        dbFieldChanged(fieldName, new Double(value), false);
	}

	public void dbFieldChanged(String fieldName, float value)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Float(value), false);
	}

	public void dbFieldChanged(String fieldName, long value)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Long(value), false);
	}

	public void dbFieldChanged(String fieldName, short value)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Short(value), false);
	}

	public void dbFieldChanged(String fieldName, Object value)
	{
		dbFieldChanged(fieldName, value, false);
	}

	public void dbFieldChanged(String fieldName, boolean value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Boolean(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, int value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Integer(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, byte value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Byte(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, char value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Character(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, float value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Float(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, double value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Double(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, long value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Long(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, short value, boolean reinitialize)
	{
		if (currentModel == null)
		{
			return;
		}
		dbFieldChanged(fieldName, new Short(value), reinitialize);
	}

	public void dbFieldChanged(String fieldName, Object value, boolean reinitialize)
	{
        if(!(this instanceof StsObject))
        {
            StsException.systemError(this, "dbFieldChanged", "for object " + StsToolkit.getSimpleClassname(this) + " is not an StsObject " + " fieldName: " + fieldName + " value: " + value);
            return;
        }
        StsObject stsObject = (StsObject)this;
        if(stsObject.getIndex() < 0)
        {
	    //    StsException.systemError("StsSerializable.dbFieldChanged() failed.  You can't persist field " + fieldName + " of nonpersistent object " + stsObject.getName());
            return;
        }
		StsChangeCmd cmd = new StsChangeCmd((StsObject)this, value, fieldName, reinitialize);
		currentModel.addTransactionCmd(toString() + "." + fieldName, cmd);
	}

	public void dbFieldChanged(String fieldName, Object oldValue, Object newValue)
	{
		if (this instanceof StsObject)
		{
			StsObject stsObject = (StsObject)this;
			if(stsObject.getIndex() < 0)
			{
				StsException.systemError("StsSerializable.dbFieldChanged() failed.  You can't persist field " + fieldName + " of nonpersistent object " + stsObject.getName());
				return;
			}
		}
		StsChangeCmd cmd = new StsChangeCmd((StsObject)this, oldValue, newValue, fieldName);
		currentModel.addTransactionCmd(toString() + "." + fieldName, cmd);
	}

	public void dbFieldChanged(String fieldName, int oldValue, int newValue)
	{
		dbFieldChanged(fieldName, new Integer(oldValue), new Integer(newValue));
	}

	public void fieldChanged(String fieldName, boolean value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, int value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, byte value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, char value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, float value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, double value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, long value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, short value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, boolean value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setBoolean(this, value);
			dbFieldChanged(fieldName, value, reinitialize);
		}
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(boolean) for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, int value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setInt(this, value);
			dbFieldChanged(fieldName, value, reinitialize);
		}
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(int) failed for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, byte value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setByte(this, value);
			dbFieldChanged(fieldName, value, reinitialize);
	   }
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(byte) for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, char value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setChar(this, value);
			dbFieldChanged(fieldName, value, reinitialize);
	   }
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(byte) for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, float value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setFloat(this, value);
		}
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(float) for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, double value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setDouble(this, value);
		}
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(double) for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, long value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setLong(this, value);
		}
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(long) for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, short value, boolean reinitialize)
	{
		try
		{
			if (currentModel == null)
			{
				return;
			}
			Field field = getField(fieldName);
			field.setAccessible(true);
			field.setShort(this, value);
		}
		catch (Exception e)
		{
			StsException.systemError("StsSerialize.fieldChanged(short) for class " + StsToolkit.getSimpleClassname(this) + " field " + fieldName);
		}
	}

	public void fieldChanged(String fieldName, Object value)
	{
		fieldChanged(fieldName, value, false);
	}

	public void fieldChanged(String fieldName, Object value, boolean reinitialize)
	{
	/*
		if (this.index() < 0)
		{
			return;
		}
	*/
		objectChanged(fieldName, value);
		dbFieldChanged(fieldName, value, reinitialize);
	}

    public Field getField(String fieldName)
	{
		StsDBTypeStsClass dbClass = (StsDBTypeStsClass)currentModel.getDatabase().getCurrentDBType(getClass().getName());
		return dbClass.getField(fieldName);
	}

	public void objectChanged(String fieldName, Object value)
	{
		try
		{
			Field field = getField(fieldName);
			if (field == null)
			{
                StsException.systemError(this, "objectChanged", "Object " + value.getClass() + " has no field: " + fieldName);
                return;
			}
			field.setAccessible(true);
			field.set(this, value);
		}
		catch (Exception e)
		{
			StsException.systemError(this, "objectChanged", "failed for field " + fieldName +
									 " of object " + value.getClass());
		}
    }

	public void setDisplayField(String fieldName, boolean value)
	{
		setDisplayField(fieldName, Boolean.toString(value));
	}

	public void setDisplayField(String fieldName, byte value)
	{
		setDisplayField(fieldName, Byte.toString(value));
	}

	public void setDisplayField(String fieldName, char value)
	{
		setDisplayField(fieldName, Character.toString(value));
	}

	public void setDisplayField(String fieldName, int value)
	{
		setDisplayField(fieldName, Integer.toString(value));
	}

	public void setDisplayField(String fieldName, float value)
	{
		setDisplayField(fieldName, Float.toString(value));
	}

	public void setDisplayField(String fieldName, double value)
	{
		setDisplayField(fieldName, Double.toString(value));
	}

	public void setDisplayField(String fieldName, long value)
	{
		setDisplayField(fieldName, Long.toString(value));
	}

	public void setDisplayField(String fieldName, short value)
	{
		setDisplayField(fieldName, Short.toString(value));
	}

	public void setDisplayField(String fieldName, Color value)
	{
		setDisplayField(fieldName, value.getRGB());
	}

	public void setDisplayField(String fieldName, String value)
	{
		try
		{
			String className = getClass().getName();
			String beanKey = className + "." + fieldName;
			currentModel.setProperty(beanKey, value);
		}
		catch(Exception e)
		{
			StsException.systemError("StsClass.setBooleanDisplayField() failed for field " + fieldName);
		}
	}

    static public String toString(Object object)
    {
        if(object == null) return "null";
        return object.toString();
    }

    static public String getLabel(StsObject object)
    {
        if(object == null) return "";
        return object.getLabel();
    }
}
