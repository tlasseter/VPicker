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

public class StsDBTypeCmd extends StsDBCommand
{
	transient private StsDBTypeClass dbTypeClass;

	public StsDBTypeCmd()
	{
		super();
	}

	public StsDBTypeCmd(StsDBTypeClass dbTypeClass)
	{
		super();
		this.dbTypeClass = dbTypeClass;
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeObject dbTypeClassType = (StsDBTypeObject)dbOutputStream.getOutputDBType(dbTypeClass);
		dbOutputStream.writeInt(dbTypeClassType.getIndex());
		dbOutputStream.writeObject(dbTypeClass, dbTypeClassType);

		if (debug)
		{
            String listName = dbOutputStream.getOutputTypeList().getName();
            System.out.println("Written definition for " + dbTypeClass.getTypeName() + " index " + dbTypeClass.getIndex() + " from outputTypeList " + listName);
		}
		dbTypeClass.flagAsWrittenToDisk();
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        int dbTypeIndex = dbInputStream.readInt();
		StsDBTypeObject dbTypeClassType = (StsDBTypeObject)dbInputStream.getInputDBType(dbTypeIndex);
		if (dbTypeClassType == null)
		{
			throw new RuntimeException("StsDBTypeCmd::read(StsDBInputStream) StsDBTypeObject is null, index = " + dbTypeIndex);
		}
		dbTypeClass = (StsDBTypeClass)dbInputStream.readObject(dbTypeClassType);

		if (debug)
		{
            String listName = dbInputStream.getInputTypeList().getName();
            System.out.println("Read definition for " + dbTypeClass.getTypeName() + " index " + dbTypeClass.getIndex() + " into inputTypeList " + listName);
		}
		dbInputStream.addStsDBTypeDefinition(dbTypeClass);
        if(debug) debugMessageRead();
    }

	public String toDebugString()
	{
		return toDebugString(dbTypeClass);
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.DB_TYPE_CMD_INDEX;
	}
}
