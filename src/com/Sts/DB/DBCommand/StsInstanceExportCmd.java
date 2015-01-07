package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsInstanceExportCmd extends StsDBCommand
{
	transient private StsObject stsObject;

	public StsInstanceExportCmd()
	{
		super();
	}

	public StsInstanceExportCmd(StsObject stsObject)
	{
		super();
		this.stsObject = stsObject;
	}

	public void abort() throws StsException
	{
		if (stsObject.getIndex() == -1) return;
		StsObjectList list = stsObject.getCreateStsClass().getStsObjectList();
		list.delete(stsObject);
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbOutputStream.getOutputDBType(stsObject);
		dbOutputStream.writeInt(dbTypeClass.getIndex());
		dbOutputStream.writeInt(stsObject.getIndex());
		dbOutputStream.writeObjectFully(stsObject, dbTypeClass);
		if (debug || true)
		{
			System.out.println("Exporting instance for " + stsObject.getClassname() + "[" + stsObject.getIndex() + "]");
		}
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
		int dbTypeIndex = dbInputStream.readInt();
		StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbInputStream.getInputDBType(dbTypeIndex);
		int objIndex = dbInputStream.readInt();
		stsObject = dbInputStream.getModelObject(dbTypeClass.getClassType(), objIndex);
		stsObject = (StsObject)dbInputStream.readObjectFully(dbTypeClass, stsObject);
		if (debug || true)
		{
			System.out.println("Importing instance for " + stsObject.getClassname() + "[" + stsObject.getIndex() + "]");
		}
        dbInputStream.addToObjects(stsObject);
	}

	public String toDebugString()
	{
		return toDebugString(stsObject);
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.INSTANCE_EXPORT_CMD_INDEX;
	}
}
