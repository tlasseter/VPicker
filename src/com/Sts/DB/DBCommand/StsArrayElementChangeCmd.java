package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.lang.reflect.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version 1.0
 */

public class StsArrayElementChangeCmd extends StsDBCommand
{
	private StsObject obj = null;
	private String fieldName;
	private Object newElement;
	private int[] arrayIndex;

	// These fields are only used in the case of the third constructor. Use with care! The change command object assumes
	// responsibility for setting and aborting the field change.
	private boolean allowRollback = false;
	private Object oldElement;

	public StsArrayElementChangeCmd()
	{
	}

	public StsArrayElementChangeCmd(StsObject obj, Object oldElement, Object newElement, String fieldName, int[] arrayIndex)
	{
		super();
		this.obj = obj;
		this.oldElement = oldElement;
		this.newElement = newElement;
		this.fieldName = fieldName;
        checkFieldName(obj, fieldName);
        this.arrayIndex = arrayIndex;
	}

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
		StsDBTypeObject indexesObjDBType = (StsDBTypeObject)dbOutputStream.getOutputDBType(int[].class);
		dbOutputStream.writeObject(arrayIndex, indexesObjDBType);
		if (newElement == null)
		{
			dbOutputStream.writeInt(StsDBTypeClass.NULL_REFERENCE);
		}
		else
		{
			StsDBTypeObject fieldObjDBType = (StsDBTypeObject)dbOutputStream.getOutputDBType(newElement);
			dbOutputStream.writeInt(fieldObjDBType.getIndex());
			dbOutputStream.writeObject(newElement, fieldObjDBType);
		}
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
		try
		{
            int index = dbInputStream.readInt();
			StsDBTypeClass objDBType = (StsDBTypeClass)dbInputStream.getInputDBType(index);
			if (objDBType == null)
			{
				throw new RuntimeException("StsChangeCmd::read(StsDBInputStream) StsDBTypeObject is null, index = " + index);
			}
			obj = (StsObject)dbInputStream.readObject(objDBType);
			int fieldNumber = dbInputStream.readInt();
			StsDBTypeObject indexesObjDBType = (StsDBTypeObject)dbInputStream.getInputDBType(int[].class);
			arrayIndex = (int[])dbInputStream.readObject(indexesObjDBType);
			index = dbInputStream.readInt();
			if (index == StsDBTypeClass.NULL_REFERENCE)
			{
				newElement = null;
			}
			else
			{
				StsDBTypeObject fieldObjDBType = (StsDBTypeObject)dbInputStream.getInputDBType(index);
				newElement = dbInputStream.readObject(fieldObjDBType);
			}
			if (obj == null)
				return;
			Field field = objDBType.getField(fieldNumber);
			if (field == null)
				return;
			field.setAccessible(true);
			Object array = field.get(obj);
			setArrayElement(array, newElement);
            if(debug) debugMessageRead();
        }
		catch (Exception e)
		{
			StsException.outputException("StsArrayElementChangeCmd::read() failed. FieldName: " + fieldName,
												  e, StsException.WARNING);
		}
	}

	private void setField(Object element)
	{
		try
		{
			StsDBTypeStsClass dbClass = (StsDBTypeStsClass)StsObject.getCurrentModel().getDatabase().getCurrentDBType(obj);
			Field field = dbClass.getField(fieldName);
			if (field == null)
			{
				return;
			}
			field.setAccessible(true);
			Object fieldObj = field.get(obj);
			if (fieldObj == null)
				return;
			setArrayElement(fieldObj, element);
		}
		catch (Exception e)
		{
			StsException.outputException("StsArrayElementChangeCmd::abort() failed. FieldName: " + fieldName,
												  e, StsException.WARNING);
		}
	}

	private void setArrayElement(Object array, Object element)
	{
		int j = arrayIndex.length - 1;
		for (int i = 0; i < j; i++)
		{
			array = Array.get(array, arrayIndex[i]);
		}
		Array.set(array, arrayIndex[j], element);
	}

	public String toDebugString()
	{
		return toDebugString(obj) + "." + fieldName + "[" + arrayIndex + "]" + toDebugString(newElement);
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.ARRAY_ELEMENT_CHANGE_CMD_INDEX;
	}
}
