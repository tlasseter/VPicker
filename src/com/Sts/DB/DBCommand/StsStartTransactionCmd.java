package com.Sts.DB.DBCommand;

import com.Sts.DB.*;

import java.io.*;
import java.util.*;

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

public class StsStartTransactionCmd extends StsDBCommand
{
	transient private String name;
	transient private long time = System.currentTimeMillis();

    public StsStartTransactionCmd()
	{
		super();
	}

	public StsStartTransactionCmd(String name)
	{
		super();
		this.name = name;
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
        StsDBFile dbFile = dbOutputStream.getDBFile();
        dbFile.setStartTransactionPosition();
        writeCmdClassIndex(dbOutputStream);
        dbOutputStream.writeUTF(name);
		dbOutputStream.writeLong(time);
        if(debug) System.out.println("write time: " + time);
        dbOutputStream.getDBFile().setLastTransactionTime(time);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
        StsDBFile dbFile = dbInputStream.getDBFile();
        dbFile.setStartTransactionPosition();
        name = dbInputStream.readUTF();
        time = dbInputStream.readLong();
		dbFile.setLastTransactionTime(time);
        if(debug)
        {
            dbFile.debugStartReadTransaction(name);
            debugMessageRead();
        }
        dbFile.incrementNTransaction();
    }

	public String toDebugString()
	{
		return "TRANSACTION " + name + "   Original Date and time: " + (new Date(time)).toString();
	}

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.START_TRANSACTION_CMD_INDEX;
	}
}