package com.Sts.DB.DBCommand;

import com.Sts.DB.*;
import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

public class StsDBMethodCmd extends StsDBCommand
{
	private Object obj = null;
	private String methodName;
	private Object[] args;

	public StsDBMethodCmd()
	{
	}

	public StsDBMethodCmd(Object obj, String methodName, Object[] args)
	{
		super();
		this.obj = obj;
		this.methodName = methodName;
		this.args = args;
	}

	public StsDBMethodCmd(Object obj, String methodName)
	{
		super();
		this.obj = obj;
		this.methodName = methodName;
		this.args = new Object[0];
	}

	public void write(StsDBOutputStream dbOutputStream) throws IOException
	{
		if (debug) debugMessageWrite();
        writeCmdClassIndex(dbOutputStream);
        StsDBTypeObject objDBType = (StsDBTypeObject)dbOutputStream.getOutputDBType(obj);
		dbOutputStream.writeInt(objDBType.getIndex());
		dbOutputStream.writeUTF(methodName);
        dbOutputStream.writeObject(obj, objDBType);
		StsDBTypeClass argsDBType = (StsDBTypeClass)dbOutputStream.getOutputDBType(args);
		dbOutputStream.writeObject(args, argsDBType);
	}

	public void read(StsDBInputStream dbInputStream) throws IOException
	{
		int index = -1;
		StsDBTypeClass argsDBType = null;
		try
		{
            index = dbInputStream.readInt();
			StsDBTypeObject objDBType = (StsDBTypeObject)dbInputStream.getInputDBType(index);
			methodName = dbInputStream.readUTF();
            obj = dbInputStream.readObject(objDBType);
			argsDBType = (StsDBTypeClass)dbInputStream.getInputDBType(Object[].class);
			args = (Object[])dbInputStream.readObject(argsDBType);
            if(debug) debugMessageRead();
            StsMethod.constructInvoke(obj, methodName, args);
		}
		catch (Exception e)
		{
            String argClassNames = getArgClassNames(args);
            StsException.outputException("StsDBMethodCmd.read() failed. Method " + methodName + " on " + obj + " with args " + argClassNames, e, StsException.WARNING);
		}
	}

    private String getArgClassNames(Object[] args)
    {
        if(args == null) return "none";
        StringBuffer stringBuffer = new StringBuffer();
        for(int n = 0; n < args.length; n++)
            stringBuffer.append(args[n].getClass().getName());
        return stringBuffer.toString();
    }

    public String toDebugString()
	{
        int nArgs = args == null ? 0 : args.length;
        return "object: " + obj + " method: " + methodName + " nArgs: " + nArgs; 
    }

	public byte getDBCommandClassIndex()
	{
		return StsDBCommand.DB_METHOD_CMD_INDEX;
	}
}
