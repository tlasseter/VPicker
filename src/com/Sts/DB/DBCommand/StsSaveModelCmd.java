package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

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

public class StsSaveModelCmd extends StsDBCommand
{
	transient private StsModel model;

	public StsSaveModelCmd()
	{
		super();
	}

	public StsSaveModelCmd(StsModel model)
	{
		super();
		this.model = model;
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        if(debug) debugMessageWrite();
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbOutputStream.getOutputDBType(model);
		dbOutputStream.writeObjectFully(model, dbTypeClass);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
		StsDBTypeObject dbTypeClass = (StsDBTypeObject)dbInputStream.getInputDBType(StsModel.class);
		model = StsObject.getCurrentModel();
		model = (StsModel)dbInputStream.readObjectFully(dbTypeClass, model);
        if (debug) debugMessageRead();
    }

	public String toDebugString()
	{
		return StsModel.class.getName();

    }

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.SAVE_MODEL_CMD_INDEX;
	}
}
