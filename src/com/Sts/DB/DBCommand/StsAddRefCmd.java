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

public class StsAddRefCmd extends StsDBCommand
{
	transient private StsObjectRefList list; // list to which ref added
	transient private StsObject obj = null;
	transient private int positionIndex;

	public StsAddRefCmd()
	{
	}

	public StsAddRefCmd(StsObjectRefList list, StsObject obj, int positionIndex)
	{
		super();
		this.list = list;
		this.obj = obj;
		this.positionIndex = positionIndex;
	}

	public void abort() throws StsException
	{
		list.deleteElement(obj);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        StsModel model = StsObject.getCurrentModel();
		int dbTypeIndex = dbInputStream.readInt();
		int objIndex = dbInputStream.readInt();
		StsDBTypeObject objDBType = (StsDBTypeObject)dbInputStream.getInputDBType(dbTypeIndex);
		obj = dbInputStream.getModelObject(objDBType.getClassType(), objIndex);
		positionIndex = dbInputStream.readInt();
		int listIndex = dbInputStream.readInt();
        if(obj == null)
            return;
        list = (StsObjectRefList)model.getStsClassObjectWithIndex(StsObjectRefList.class, listIndex);
        if(debug) debugMessageRead();
        if (list == null)
		{
			StsException.systemError(this, "read", "Could not find StsObjectRefList with index " + listIndex);
			return;
		}
		if (!list.insertObject(obj, positionIndex))
		{
			StsException.systemError(this, "read", toString()  + " Failed inserting object.");
			return;
		}
    }

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
		if (debug) debugMessageWrite();
        if(list.getIndex() == -1 || obj.getIndex() == -1)
            return;        
        writeCmdClassIndex(dbOutputStream);
        StsDBType dbType = dbOutputStream.getOutputDBType(obj);
		dbOutputStream.writeInt(dbType.getIndex());
		dbOutputStream.writeInt(obj.getIndex());
		dbOutputStream.writeInt(positionIndex);
		dbOutputStream.writeInt(list.getIndex());
	}

	public String toDebugString()
	{
        String listString = list == null ? "null" : list.toString();
        return "StsObjectRefList " + listString + " added Object " + obj.toDebugString() + " in position " + positionIndex;
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.ADD_REF_CMD_INDEX;
	}
}
