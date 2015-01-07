package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
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

public class StsInstanceDeleteCmd extends StsDBCommand
{
	transient private StsObject obj = null;
    transient private int objectIndex = -99;
    transient private Class c;

	public StsInstanceDeleteCmd()
	{
		super();
	}

	public StsInstanceDeleteCmd(StsObject obj)
	{
		super();
		this.obj = obj;
        objectIndex = obj.getIndex();
        c = obj.getClass();
	}

	public void abort() throws StsException
	{
		StsObject stsObject = (StsObject)obj;
		stsObject.getCreateStsClass().getStsObjectList().add(obj);
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        if(objectIndex == -1) return; //  && obj.getClass() != StsProject.class) return;
        writeCmdClassIndex(dbOutputStream);
		StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbOutputStream.getOutputDBType(obj);
		dbOutputStream.writeInt(dbTypeClass.getIndex());
		dbOutputStream.writeInt(objectIndex);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        int dbTypeIndex = dbInputStream.readInt();
		int objectIndex = dbInputStream.readInt();
		StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbInputStream.getInputDBType(dbTypeIndex);
		obj = dbInputStream.getModelObject(dbTypeClass.getClassType(), objectIndex);
		c = obj.getClass();
		StsModel model = StsObject.getCurrentModel();
		model.deleteStsClassObject(c, objectIndex);
        dbInputStream.deleteObjectFromList(obj);
        if(debug) debugMessageRead();
    }

	public String toDebugString()
	{
		return toDebugString(c, obj.getIndex());
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.INSTANCE_DELETE_CMD_INDEX;
	}
}
