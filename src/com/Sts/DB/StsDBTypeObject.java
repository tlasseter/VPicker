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

public abstract class StsDBTypeObject extends StsDBType implements StsSerializable
{
	transient private StsDBObjectIO objectIO;

	public StsDBTypeObject()
	{
	}

	public StsDBTypeObject(int index, Class classType)
	{
		super(index, classType);
		this.objectIO = getObjectIO();
	}

	public StsDBTypeObject(int index)
	{
		super(index);
		this.objectIO = getObjectIO();
	}

    public StsDBTypeObject(int index, Class classType, StsDBObjectIO objectIO)
	{
		super(index, classType);
		this.objectIO = objectIO;
	}

	public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
		objectIO.exportObject(objectTrader, obj);
	}

	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		objectIO.writeObject(out, obj);
	}

	public void writeObjectFully(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		objectIO.writeObject(out, obj);
	}

	public Object readObject(StsDBInputStream in, Object object) throws IOException
	{
		return objectIO.readObject(in, object);
	}

	public Object readObject(StsDBInputStream in) throws IOException
	{
		return objectIO.readObject(in);
	}

	public Object readObjectFully(StsDBInputStream in, Object object) throws IOException
	{
		return objectIO.readObject(in, object);
	}

	public Object readObjectFully(StsDBInputStream in) throws IOException
	{
		return objectIO.readObject(in);
	}

	public Object copyObject(Object oldObject) throws IOException, IllegalAccessException
	{
		return objectIO.copyObject(oldObject);
	}

	protected void initializeAfterLoad(StsDBObjectTypeList typeList) throws ClassNotFoundException
	{
		super.initializeAfterLoad(typeList);
		this.objectIO = getObjectIO();
	}

	protected void initializeObjectIO()
	{
		this.objectIO = getObjectIO();
	}

	protected StsDBObjectIO getObjectIO()
	{
		if (StsCustomSerializable.class.isAssignableFrom(this.getClassType()))
		{
			return new StsDBObjectIOCustom(this);
		}
		return new StsDBObjectIOObject(this);
	}

	protected abstract Object newInstance();
}
