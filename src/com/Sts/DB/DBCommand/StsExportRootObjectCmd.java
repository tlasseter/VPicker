package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;

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

public class StsExportRootObjectCmd extends StsDBCommand
{
	static private StsObject rootObject;

	public StsExportRootObjectCmd()
	{
		super();
	}

	public StsExportRootObjectCmd(StsObject rootObject)
	{
		super();
		this.rootObject = rootObject;
	}

	static public StsObject getRootObject()
	{
		return rootObject;
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbOutputStream.getOutputDBType(rootObject);
		dbOutputStream.writeInt(dbTypeClass.getIndex());
		dbOutputStream.writeInt(rootObject.getIndex());
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        int dbTypeIndex = dbInputStream.readInt();
		StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbInputStream.getInputDBType(dbTypeIndex);
		int objIndex = dbInputStream.readInt();

		rootObject = dbInputStream.getModelObject(dbTypeClass.getClassType(), objIndex);
        if(debug) debugMessageRead();
	}

	public String toDebugString()
	{
		return toDebugString(rootObject);
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.EXPORT_ROOT_OBJECT_CMD_INDEX;
	}
}
