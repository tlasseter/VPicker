package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.lang.reflect.*;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

public class StsArrayChangeCmd extends StsDBCommand
{
	private StsObject obj = null;
	private String fieldName;
	private Object element;
	private int arrayIndex;
	private boolean reinitialize = false;

	// These fields are only used in the case of the third constructor. Use with care! The change command object assumes
	// responsibility for setting and aborting the field change.
	private boolean allowRollback = false;
	private Object oldElement;

	public StsArrayChangeCmd()
	{
	}

	public StsArrayChangeCmd(StsObject obj, Object element, String fieldName, int arrayIndex, boolean reinitialize)
	{
		super();
		this.obj = obj;
		this.element = element;
		this.fieldName = fieldName;
        checkFieldName(obj, fieldName);
        this.arrayIndex = arrayIndex;
        if(arrayIndex < 0)
            System.out.println("StsArrayChangeCmd error. bad index: " + arrayIndex);
        this.reinitialize = reinitialize;
	}
/*
	public StsArrayChangeCmd(StsObject obj, Object element, String fieldName)
	{
		super();
		this.obj = obj;
		this.element = element;
		this.fieldName = fieldName;
		reinitialize = false;
	}

	public StsArrayChangeCmd(StsObject obj, Object oldElement, Object element, String fieldName)
	{
		super();
		this.obj = obj;
		this.oldElement = oldElement;
		this.element = element;
		this.fieldName = fieldName;
		allowRollback = true;
		reinitialize = false;
	}
*/
	public void abort() throws StsException
	{
		if (allowRollback)
			setField(oldElement);
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        if(!obj.isPersistent()) return;
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeClass objDBType = (StsDBTypeClass)dbOutputStream.getOutputDBType(obj);
		dbOutputStream.writeInt(objDBType.getIndex());
		dbOutputStream.writeObject(obj, objDBType);
		int fieldNumber = objDBType.getIndexOfField(fieldName);
		dbOutputStream.writeInt(fieldNumber);
		if (element == null)
		{
			dbOutputStream.writeInt(StsDBTypeClass.NULL_REFERENCE);
		}
		else
		{
			StsDBTypeObject fieldObjDBType = (StsDBTypeObject)dbOutputStream.getOutputDBType(element);
			dbOutputStream.writeInt(fieldObjDBType.getIndex());
			dbOutputStream.writeInt(arrayIndex);
			dbOutputStream.writeObject(element, fieldObjDBType);
		}
		dbOutputStream.writeBoolean(reinitialize);
	}

	private void setField(Object element)
	{
		try
		{
			StsDBTypeStsClass dbClass = (StsDBTypeStsClass)StsObject.getCurrentModel().getDatabase().getCurrentDBType(obj);
			//Field field = obj.getClass().getField(fieldName);
			Field field = dbClass.getField(fieldName);
			if (field == null)
			{
				return;
			}
			field.setAccessible(true);
			Object fieldObj = field.get(obj);
			if (fieldObj == null)
				return;
			Object[] array = (Object[])fieldObj;
			array[arrayIndex] = element;
		}
		catch (Exception e)
		{
			StsException.outputException("StsChangeCmd::abort() failed. FieldName: " + fieldName,
												  e, StsException.WARNING);
		}
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        int index;
        StsDBTypeClass objDBType;
        StsDBTypeObject fieldObjDBType;
        Field field;
        Object fieldObj;
        Object[] array;
        
        try
		{
            index = dbInputStream.readInt();
			objDBType = (StsDBTypeClass)dbInputStream.getInputDBType(index);
			if (objDBType == null)
			{
				throw new RuntimeException("StsArrayChangeCmd::read(StsDBInputStream) StsDBTypeObject is null, index = " + index);
			}
			obj = (StsObject)dbInputStream.readObject(objDBType);
			int fieldNumber = dbInputStream.readInt();
			index = dbInputStream.readInt();
			if (index == StsDBTypeClass.NULL_REFERENCE)
			{
				element = null;
			}
			else
			{
				arrayIndex = dbInputStream.readInt();
				fieldObjDBType = (StsDBTypeObject)dbInputStream.getInputDBType(index);
				element = dbInputStream.readObject(fieldObjDBType);
			}
			reinitialize = dbInputStream.readBoolean();
			if (obj == null)
				return;
			if (element == null)
				return;
			field = objDBType.getField(fieldNumber);
			if (field == null)
				return;
			field.setAccessible(true);
			fieldObj = field.get(obj);
			if (fieldObj == null)
				return;
			array = (Object[])fieldObj;
			if (array.length < arrayIndex + 1)
            {
                array = (Object[])StsMath.arrayAddElement(array, element);
                field.set(obj, array);
            }
            else
			    array[arrayIndex] = element;

			if (obj != null && reinitialize)
				dbInputStream.addToObjects((obj));
            if(debug) debugMessageRead();
        }
		catch (Exception e)
		{
			StsException.outputException("StsArrayChangeCmd::read() failed. FieldName: " + fieldName,
												  e, StsException.WARNING);
		}
	}

	public String toDebugString()
	{
		return toDebugString(obj) + "." + fieldName + "[" + arrayIndex + "]" + toDebugString(element);
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.ARRAY_CHANGE_CMD_INDEX;
	}
}
