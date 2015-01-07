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

public class StsDeleteRefCmd extends StsDBCommand
{
	transient private StsObjectRefList list;
    transient private int listIndex = -1;
    transient private StsObject obj = null;
    transient private int positionIndex;

	public StsDeleteRefCmd()
	{
	}

	public StsDeleteRefCmd(StsObjectRefList list, StsObject obj, int positionIndex)
	{
		super();
		this.list = list;
        listIndex = list.getIndex();
        this.obj = obj;
		this.positionIndex = positionIndex;
	}

	public void abort() throws StsException
	{
		list.add(obj);
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if (debug) debugMessageWrite();
        if(listIndex == -1) return;
        writeCmdClassIndex(dbOutputStream);
        dbOutputStream.writeInt(listIndex);
		dbOutputStream.writeInt(positionIndex);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        int listIndex = dbInputStream.readInt();
		positionIndex = dbInputStream.readInt();

		StsModel model = StsObject.getCurrentModel();
		StsObjectRefList refList = (StsObjectRefList)model.getStsClassObjectWithIndex(StsObjectRefList.class, listIndex);
		if (refList == null)
		{
			StsException.systemError(this, "read", "Could not find " + toString());
			return;
		}
		if (!refList.delete(positionIndex))
		{
			StsException.systemError(this, "read", "Couldn't delete element  " + positionIndex + " in objectRefList " + refList.getName());
			return;
		}
        if(debug) debugMessageRead();
    }

	public String toDebugString()
	{
        String listString = list == null ? "null" : list.toString();
        return "StsObjectRefList " + listString + " deleted Object in position " + positionIndex;
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.DELETE_REF_CMD_INDEX;
	}
}
