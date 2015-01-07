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

public class StsAddTransientModelObjectCmd extends StsDBCommand
{
	private String fieldName;
	private Object fieldObj;
	private boolean reinitialize = false;

	// These fields are only used in the case of the third constructor. Use with care! The change command object assumes
	// responsibility for setting and aborting the field change.
	private boolean allowRollback = false;
	private Object oldFieldObject;

	public StsAddTransientModelObjectCmd()
	{
	}

	public StsAddTransientModelObjectCmd(Object fieldObj, String fieldName, boolean reinitialize)
	{
		super();
        if(fieldName == null)
            StsException.systemError(this, "constructor", "fieldName field is null.");
        super.checkFieldName(StsModel.class, fieldName);
        this.fieldObj = fieldObj;
		this.fieldName = fieldName;
		this.reinitialize = reinitialize;
	}

	public StsAddTransientModelObjectCmd(Object fieldObj, String fieldName)
	{
		super();
        if(fieldName == null) StsException.systemError(this, "constructor", "fieldName field is null.");
		this.fieldObj = fieldObj;
		this.fieldName = fieldName;
		reinitialize = false;
	}

	public StsAddTransientModelObjectCmd(Object oldFieldObject, Object fieldObj, String fieldName)
	{
		super();
        if(fieldName == null)
            StsException.systemError(this, "constructor", "fieldName field is null.");
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
        try
        {
            if (debug) debugMessageWrite();
            writeCmdClassIndex(dbOutputStream);

            StsDBTypeClass.setTypeIsTemporary(true);
//            boolean saveDebug = dbOutputStream.debug;
//            dbOutputStream.debug = true;
            if (fieldObj == null)
            {
                dbOutputStream.writeInt(StsDBTypeClass.NULL_REFERENCE);
            }
            else
            {
                dbOutputStream.writeUTF(fieldName);
                StsDBTypeObject fieldObjDBType = (StsDBTypeObject)dbOutputStream.getOutputDBType(fieldObj);
                dbOutputStream.writeInt(fieldObjDBType.getIndex());
                dbOutputStream.writeObject(fieldObj, fieldObjDBType);
            }
//            dbOutputStream.debug = saveDebug;
            dbOutputStream.writeBoolean(reinitialize);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "write()", e);
        }
        finally
        {
            StsDBTypeClass.setTypeIsTemporary(false);
        }
    }

	private void setField(Object newFieldObject)
	{
		try
		{
            Field field = StsModel.class.getField(fieldName);
			if (field == null) return;
			field.setAccessible(true);
            StsModel model = StsModel.getCurrentModel();
            field.set(model, newFieldObject);
		}
		catch (Exception e)
		{
			StsException.outputException("StsAddTransientModelObjectCmd::abort() failed. FieldName: " + fieldName,
												  e, StsException.WARNING);
		}
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
//        StsDBFile dbFile = null;
        try
		{
            fieldName = dbInputStream.readUTF();
            Field field = StsModel.class.getField(fieldName);
            StsModel model = StsObject.getCurrentModel();

            StsDBTypeClass.setTypeIsTemporary(true);

            int index = dbInputStream.readInt();
			StsDBTypeObject fieldObjDBType;
			if (index == StsDBTypeClass.NULL_REFERENCE)
			{
				fieldObj = null;
			}
			else
			{
                fieldObj = field.get(model);
                fieldObjDBType = (StsDBTypeObject)dbInputStream.getInputDBType(index);
                fieldObj = fieldObjDBType.readObject(dbInputStream, fieldObj);
			}
			reinitialize = dbInputStream.readBoolean();
			if (model == null) return;
			if (fieldObj == null) return;
            if (field == null) return;
            if(debug) debugMessageRead();
			field.setAccessible(true);
			field.set(model, fieldObj);
		}
		catch (Exception e)
		{
			StsException.outputException("StsAddTransientModelObjectCmd::read() failed. FieldName: " + fieldName,
										 e, StsException.WARNING);
		}
        finally
        {
            StsDBTypeClass.setTypeIsTemporary(false);
        }
    }

	public String toDebugString()
	{
		return fieldName + " " + toDebugString(fieldObj);
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.ADD_TRANS_MODEL_OBJECT_CMD_INDEX;
	}
}