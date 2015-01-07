package com.Sts.DB.DBCommand;

import com.Sts.DB.*;

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

public class StsSimpleObjectAddCmd extends StsDBCommand
{
	// making this static so that it can be set before writing out and read after reading in
	// without having to have a reference to 'this' object
	private static Object staticObject;

	public StsSimpleObjectAddCmd()
	{
		super();
	}

	public StsSimpleObjectAddCmd(Object obj)
	{
		super();
		staticObject = obj;
	}

	public static void setObject(Object obj)
	{
		staticObject = obj;
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbOutputStream.getOutputDBType(staticObject);
		dbOutputStream.writeInt(dbTypeClass.getIndex());
		dbOutputStream.writeObjectFully(staticObject, dbTypeClass);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        int dbTypeIndex = dbInputStream.readInt();
		StsDBTypeObject dbTypeClassType = (StsDBTypeObject)dbInputStream.getInputDBType(dbTypeIndex);
		staticObject = dbInputStream.readObjectFully(dbTypeClassType, staticObject);
        if(debug) debugMessageRead();
    }

	public String toDebugString()
	{
		if (staticObject != null)
			return staticObject.getClass().getName();
		else
			return "null object";
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.SIMPLE_OBJECT_ADD_CMD_INDEX;
	}
}
