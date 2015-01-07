package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
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

public class StsChangeCmd extends StsDBCommand
{
	private Object obj = null;
	private String fieldName;
	private Object fieldObj;
	private boolean reinitialize = false;

	// These fields are only used in the case of the third constructor. Use with care! The change command object assumes
	// responsibility for setting and aborting the field change.
	private boolean allowRollback = false;
	private Object oldFieldObject;

	public StsChangeCmd()
	{
	}

	public StsChangeCmd(StsObject obj, Object fieldObj, String fieldName, boolean reinitialize)
	{
		super();
        if(fieldName == null)
            StsException.systemError(this, "constructor", "fieldName field is null.");
        checkFieldName(obj, fieldName);
        this.obj = obj;
		this.fieldObj = fieldObj;
		this.fieldName = fieldName;
        checkFieldName(obj, fieldName);
        this.reinitialize = reinitialize;
	}

	public StsChangeCmd(StsModel obj, Object fieldObj, String fieldName)
	{
		super();
        if(fieldName == null)
            StsException.systemError(this, "constructor", "fieldName field is null.");
        this.obj = obj;
		this.fieldObj = fieldObj;
		this.fieldName = fieldName;
		reinitialize = false;
	}

	public StsChangeCmd(StsObject obj, Object oldFieldObject, Object fieldObj, String fieldName)
	{
		super();
        if(fieldName == null)
            StsException.systemError(this, "constructor", "fieldName field is null.");
        this.obj = obj;
		this.oldFieldObject = oldFieldObject;
		this.fieldObj = fieldObj;
		this.fieldName = fieldName;
		allowRollback = true;
		setField(fieldObj);
		reinitialize = false;
	}

	public void abort() throws StsException
	{
		if (allowRollback)
		{
			setField(oldFieldObject);
		}
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        StsDBTypeClass objDBType = (StsDBTypeClass)dbOutputStream.getOutputDBType(obj);
        int fieldNumber = objDBType.getIndexOfField(fieldName);
        if(fieldNumber == -1)
        {
            String objClassname = StsToolkit.getSimpleClassname(objDBType.getClassType());
            StsException.systemError(this, "write", objClassname + "." + fieldName + " not found. Perhaps it is transient. Aborting write.");
            return;
        }
        writeCmdClassIndex(dbOutputStream);
        dbOutputStream.writeInt(objDBType.getIndex());
		dbOutputStream.writeObject(obj, objDBType);

		dbOutputStream.writeInt(fieldNumber);

        if (fieldObj == null)
		{
			dbOutputStream.writeInt(StsDBTypeClass.NULL_REFERENCE);
		}
		else
		{
			StsDBTypeObject fieldObjDBType = (StsDBTypeObject)dbOutputStream.getOutputDBType(fieldObj);
			dbOutputStream.writeInt(fieldObjDBType.getIndex());
			dbOutputStream.writeObject(fieldObj, fieldObjDBType);
		}
		dbOutputStream.writeBoolean(reinitialize);
	}

	private void setField(Object newFieldObject)
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
			field.set(obj, newFieldObject);
		}
		catch (Exception e)
		{
			StsException.outputException("StsChangeCmd::abort() failed. FieldName: " + fieldName,
												  e, StsException.WARNING);
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
			obj = dbInputStream.readObject(objDBType);
			int fieldNumber = dbInputStream.readInt();
			index = dbInputStream.readInt();
			StsDBTypeObject fieldObjDBType = null;
			if (index == StsDBTypeClass.NULL_REFERENCE)
			{
				fieldObj = null;
			}
			else
			{
				fieldObjDBType = (StsDBTypeObject)dbInputStream.getInputDBType(index);
				fieldObj = dbInputStream.readObject(fieldObjDBType);
			}
			reinitialize = dbInputStream.readBoolean();
			if (obj == null)
				return;
			if (fieldObj == null)
				return;
			Field field = objDBType.getField(fieldNumber);
			if (field == null)
				return;
			field.setAccessible(true);
			field.set(obj, fieldObj);
            if(debug) debugMessageRead();
//			if (obj != null && reinitialize)
//				dbInputStream.addToObjects(((StsObject)obj));
		}
		catch (Exception e)
		{
			StsException.outputException("StsChangeCmd::read() failed. FieldName: " + fieldName,
										 e, StsException.WARNING);
		}
	}

	public String toDebugString()
	{
		return toDebugString(obj) + "." + fieldName + " " + toDebugString(fieldObj);
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.CHANGE_CMD_INDEX;
	}
}
