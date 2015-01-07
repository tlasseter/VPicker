package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
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

public class StsEndDBCmd extends StsDBCommand
{
    static final boolean debug = false;

    public StsEndDBCmd()
	{
		super();
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
    {
		dbInputStream.setPositionEndDB();
        // long position = dbInputStream.readLong();
        // long currentPosition = StsModel.getCurrentModel().getDatabase().getPosition();
        if(Main.isDbCmdDebug) System.out.println("    EndDBCmd read");
        // if(Main.isDbDebug) System.out.println("    EndDBCmd read at position: " + position + " currentPosition: " + currentPosition);
        // dbInputStream.endDB();
    }

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
		dbOutputStream.endDB();
        writeCmdClassIndex(dbOutputStream);
        // long position = StsModel.getCurrentModel().getDatabase().getPosition();
        // dbOutputStream.writeLong(position);
        if(Main.isDbDebug) System.out.println("    EndDBCmd written");
    }

	public String toDebugString()
	{
		return "End DB command";
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.END_DB_CMD_INDEX;
	}
}
