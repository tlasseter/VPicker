package com.Sts.DB;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

import java.io.*;
import java.lang.reflect.*;

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

public class StsDBTypeModel extends StsDBTypeClass
{
	public StsDBTypeModel()
	{
		super();
	}

	public StsDBTypeModel(int index)
	{
		super(index, StsModel.class);
	}

	public Object copyArrayValues(Object arrayObject) throws IOException, IllegalAccessException
	{
		throw new RuntimeException("This method should never be called");
	}

	public void copyField(Object oldObject, Object newObject, Field field) throws IOException, IllegalAccessException
	{
		throw new RuntimeException("This method should never be called");
	}

	protected Object newInstance()
	{
		throw new RuntimeException("This method should never be called");
	}

	public Object readArrayValues(StsDBInputStream in, Object arrayObject)
	{
		throw new RuntimeException("This method should never be called");
	}

	public void readField(StsDBInputStream in, Object classObject, Field field) throws IllegalAccessException, IOException
	{
		in.readInt();
		if (field != null)
			field.set(classObject, StsObject.getCurrentModel());
	}

	public void writeArrayValues(StsDBOutputStream out, Object arrayObject) throws IOException, IllegalAccessException
	{
		throw new RuntimeException("This method should never be called");
	}

	public void writeField(StsDBOutputStream out, Object classObject, Field field) throws IllegalAccessException, IOException
	{
		out.writeInt(getIndex());
	}

	public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
	}

	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
	}

	public Object readObject(StsDBInputStream in, Object object) throws IOException
	{
		return StsObject.getCurrentModel();
	}

	public Object readObject(StsDBInputStream in) throws IOException
	{
		return StsObject.getCurrentModel();
	}

	public Object copyObject(Object oldObject) throws IOException, IllegalAccessException
	{
		return StsObject.getCurrentModel();
	}
}
