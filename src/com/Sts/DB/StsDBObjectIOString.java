package com.Sts.DB;

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

public class StsDBObjectIOString implements StsDBObjectIO
{
	public StsDBObjectIOString()
	{
	}

	public void writeObject(StsDBOutputStream out, Object obj) throws IOException
	{
		out.writeBoolean(obj == null);
		if (obj == null)
		{
			return;
		}
		out.writeUTF((String)obj);
	}

	public Object readObject(StsDBInputStream in) throws IOException
	{
		boolean isNull = in.readBoolean();
		if (isNull)
		{
			return null;
		}
		return in.readUTF();
	}

	/**
	 * readObject
	 *
	 * @param in StsDBInputStream
	 * @param obj Object
	 * @return Object
	 * @throws IOException
	 */
	public Object readObject(StsDBInputStream in, Object obj) throws IOException
	{
		return readObject(in);
	}

	public Object copyObject(Object oldObject) throws IllegalAccessException
	{
		return new String((String)oldObject);
	}

	public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
	}
}
