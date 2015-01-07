package com.Sts.DB;

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

public class StsDBOutputStream extends ObjectOutputStream
{
	public static boolean debug = Main.isDbCmdDebug;

	private StsDBObjectTypeList outputDBClasses = null;
	private StsDBFile dbFile = null;

	public StsDBOutputStream(StsDBFile dbFile, OutputStream out, StsDBObjectTypeList outputDBClasses) throws IOException
	{
		super(out);
		this.dbFile = dbFile;
		this.outputDBClasses = outputDBClasses;
	}

	public StsDBFile getDBFile()
	{
		return dbFile;
	}

	public final void writeObject(Object obj, StsDBTypeObject dbClass)
	{
		try
		{
			dbClass.writeObject(this, obj);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("Exception in StsDBOutputStream.write()\n" + e);
		}
	}

	public final void writeObjectFully(Object obj, StsDBTypeObject dbClass)
	{
		try
		{
			dbClass.writeObjectFully(this, obj);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("Exception in StsDBOutputStream.write()\n" + e);
		}
	}

	public StsDBObjectTypeList getOutputTypeList()
	{
		return outputDBClasses;
	}

	public StsDBType getOutputDBType(Class c)
	{
		return outputDBClasses.getDBType(c);
	}

	public StsDBType getOutputDBType(Object o)
	{
		return outputDBClasses.getDBType(o);
	}
/*
	public Object[] getNonDefaultDBClassArray()
	{
		return outputDBClasses.getNonDefaultDBClassArray();
	}
*/
	public void writeBoolean(boolean value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeBoolean(value);
	}

	public void writeByte(int value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeByte(value);
	}

	public void writeChar(char value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeChar(value);
	}

	public void writeDouble(double value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeDouble(value);
	}

	public void writeFloat(float value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeFloat(value);
	}

	public void writeInt(int value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeInt(value);
	}

	public void writeLong(long value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeLong(value);
	}

	public void writeShort(short value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeShort(value);
	}

	public void writeUTF(String value) throws java.io.IOException
	{
		if (debug) System.out.println("Output---->" + value);
		super.writeUTF(value);
	}

    public void endDB()
    {
        dbFile.setPositionEndDB();
    }
	public static void main(String[] args)
	{
		try
		{
			String filename = "c:\\temp\\test.db";
			FileOutputStream fos = new FileOutputStream(filename);
			StsDBOutputStream out = new StsDBOutputStream(null, fos, null);
/*
			String s = "test";
			out.write(s);
			StsColor color = new StsColor(1.0f, 1.0f, 1.0f);
			out.write(color);
*/
			out.flush();
			fos.close();
		}
		catch(Exception exc)
		{
			System.out.println("Exception in DBOutputStream.main()\n" + exc);
		}
	}
}

