package com.Sts.DB;

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

public class StsDBObjectIOCustom extends StsDBObjectIOObject
{
	public StsDBObjectIOCustom(StsDBTypeObject classType)
	{
		super(classType);
	}

	/** copyObject
	  * @return Object
	  */
	 public Object copyObject(Object obj) throws IllegalAccessException, IOException
	 {
		 Object newObject = super.copyObject(obj);
		 return newObject;
	 }

	 /**
	  * readObject
	  *
	  * @param in StsDBInputStream
	  * @return Object
	  * @throws IOException
	  */
	 public Object readObject(StsDBInputStream in) throws IOException
	 {
		 return this.readObject(in, null);
	 }

	/**
	 * readObject
	 *
	 * @param in StsDBInputStream
	 * @param object Object
	 * @return Object
	 * @throws IOException
	 */
	public Object readObject(StsDBInputStream in, Object object) throws IOException
	{
		try
		{
			Object obj = super.readObject(in, object);
			((StsCustomSerializable)obj).readObject(in);
			return obj;
		}
		catch (Exception e)
		{
			StsException.outputException("StsDBObjectIOCustom.read(in, obj) failed.", e, StsException.WARNING);
			return null;
		}
	}

	/**
	 * writeObject
	 *
	 * @param out StsDBOutputStream
	 * @param obj Object
	 * @throws IOException
	 * @throws IllegalAccessException
	 */
	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		super.writeObject(out, obj);
		((StsCustomSerializable)obj).writeObject(out);
	}


    public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
    {
		super.exportObject(objectTrader, obj);
		((StsCustomSerializable)obj).exportObject(objectTrader);
    }
}
