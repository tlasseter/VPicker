package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: this command adds a new object or completely replaces an existing object</p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsInstanceAddCmd extends StsDBCommand
{
	transient private StsObject stsObject = null;

	public StsInstanceAddCmd()
	{
		super();
	}

	public StsInstanceAddCmd(StsObject obj)
	{
		super();
		this.stsObject = obj;
	}

	public void abort() throws StsException
	{
		if (stsObject.getIndex() == -1)
			return;
		StsObjectList list = stsObject.getCreateStsClass().getStsObjectList();
		list.delete(stsObject);
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        int index = stsObject.getIndex();
        if(index == -1) return; //  && stsObject.getClass() != StsProject.class) return;
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbOutputStream.getOutputDBType(stsObject);
		dbOutputStream.writeInt(dbTypeClass.getIndex());
		dbOutputStream.writeInt(stsObject.getIndex());
		dbOutputStream.writeObjectFully(stsObject, dbTypeClass);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        int dbTypeIndex = dbInputStream.readInt();
		StsDBTypeClass dbTypeClass = (StsDBTypeClass)dbInputStream.getInputDBType(dbTypeIndex);
		int objIndex = dbInputStream.readInt();

		Object existingObject = dbInputStream.getModelObjectOrNull(dbTypeClass.getClassType(), objIndex);
		if (existingObject == null)
		{
			stsObject = dbInputStream.getModelObject(dbTypeClass.getClassType(), objIndex);
			if (StsDBType.debug)
			{
				dbTypeClass.sanityCheck(stsObject);
			}
			stsObject = (StsObject)dbInputStream.readObjectFully(dbTypeClass, stsObject);
		}
		else
		{
			stsObject = (StsObject)dbInputStream.readObjectFully(dbTypeClass, existingObject);
		}
        if(debug) debugMessageRead();
		dbInputStream.addToObjects(stsObject);
	}

	public String toDebugString()
	{
		return stsObject.toDebugString();
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.INSTANCE_ADD_CMD_INDEX;
	}
}
