package com.Sts.DB;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

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

public class StsDBInputStream extends ObjectInputStream
{
	static final public boolean debug = Main.isDbIODebug;
	private InputDBClasses inputDBClasses = null;
	private ArrayList objects = null;
	private StsDBObjectTypeList currentDBClasses = null;
	private StsDBFile dbFile = null;
    public boolean endDBFound = false;

	public StsDBInputStream(StsDBFile dbFile, InputStream in) throws IOException
	{
		super(in);
		this.dbFile = dbFile;
		inputDBClasses = new InputDBClasses();
		currentDBClasses = new StsDBObjectTypeList("inputDBClasses");
		objects = new ArrayList(1000);
	}

	public StsDBInputStream(StsDBFile dbFile, InputStream in, StsDBInputStream oldDBInputStream) throws IOException
	{
		super(in);
		this.dbFile = dbFile;
		if (oldDBInputStream == null || oldDBInputStream.objects == null)
		{
			objects = new ArrayList(1000);
		}
		else
		{
			objects = oldDBInputStream.objects;
		}
		if (oldDBInputStream == null || oldDBInputStream.currentDBClasses == null)
		{
			currentDBClasses = new StsDBObjectTypeList("inputDBClasses");
		}
		else
		{
			currentDBClasses = oldDBInputStream.currentDBClasses;
		}
		if (oldDBInputStream == null || oldDBInputStream.inputDBClasses == null)
		{
			inputDBClasses = new InputDBClasses();
		}
		else
		{
			inputDBClasses = oldDBInputStream.inputDBClasses;
		}
	}

	public StsDBInputStream(StsDBFile dbFile, InputStream in, StsDBInputStream oldDBInputStream, StsDBObjectTypeList currentDBClasses) throws IOException
	{
		super(in);
		this.dbFile = dbFile;
		if (oldDBInputStream == null || oldDBInputStream.objects == null)
		{
			objects = new ArrayList(1000);
		}
		else
		{
			objects = oldDBInputStream.objects;
		}
		this.currentDBClasses = currentDBClasses;
		if (oldDBInputStream == null || oldDBInputStream.inputDBClasses == null)
		{
			inputDBClasses = new InputDBClasses();
		}
		else
		{
			inputDBClasses = oldDBInputStream.inputDBClasses;
		}
	}

	public StsDBFile getDBFile()
	{
		return dbFile;
	}

	public void addToObjects(StsObject object)
	{
		if (object != null)
			objects.add(object);
	}

    public void deleteObjectFromList(StsObject deleteObject)
    {
        Iterator objectIter = objects.iterator();
        while(objectIter.hasNext())
        {
            Object object = objectIter.next();
            if(object == deleteObject)
                objectIter.remove();
        }
    }

    public void addStsDBTypeDefinition(StsDBTypeClass dbTypeClass)
	{
		try
		{
			dbTypeClass.initializeAfterLoad(inputDBClasses);
			StsDBTypeClass currentDBClass = (StsDBTypeClass)currentDBClasses.getDBTypeForceInsert(dbTypeClass);
			inputDBClasses.setDBClass(currentDBClass, dbTypeClass);
		}
		catch (ClassNotFoundException e)
		{
			StsException.outputException("StsDBFile.addStsDBTypeDefintion() failed.", e, StsException.WARNING);
		}
	}

	protected void initializeObjects(StsModel model)
	{
		if (objects == null)
			return;
		int nObjects = objects.size();
		if (nObjects == 0)
			return;

		boolean initialized = false;
		int iter = 0;
		int maxIter = 3;
		// model.project.setDatabaseInfo(model.getDatabase());
 //       model.project.initialize(model);
        while (!initialized && iter++ < maxIter)
		{
			initialized = true;
			for (int n = 0; n < nObjects; n++)
			{
				StsObject obj = (StsObject)objects.get(n);
                if(!obj.isPersistent()) continue;
                if (!obj.initialize(model))
				{
					// obj.delete();
					initialized = false;
					if (iter == maxIter)
//                    if(mainDebug && iter == maxIter)
						System.out.println("     failed to classInitialize " + obj.getName());
				}
				else if (debug)
					System.out.println("Initialized " + obj.getName());
			}
		}
		objects.clear();
	}

	public StsObject getModelObject(Class c, int index)
	{
		StsModel model = StsObject.getCurrentModel();
		if (index < 0)
			return null;
		else
		{
			StsObject object;
			object = model.getStsClassObjectWithIndex(c, index);
			if (object != null)
				return object;
			return model.getEmptyStsClassObject(c, index);
		}
	}

	public StsObject getModelObjectOrNull(Class c, int index)
	{
		StsModel model = StsObject.getCurrentModel();
		if (index < 0)
			return null;
		else
		{
			StsObject object;
			object = model.getStsClassObjectWithIndex(c, index);
			return object;
		}
	}

	public Object readObject(StsDBTypeObject dbClass)
	{
		try
		{
			return dbClass.readObject(this);
		}
		catch (Exception e)
		{
			StsException.outputException("DBInputStream::readObject(Class, Object) failed for object of class " + dbClass,
				e, StsException.WARNING);
			return null;
		}
	}

	public Object readObject(StsDBTypeObject dbClass, Object object)
	{
		try
		{
			return dbClass.readObject(this, object);
		}
		catch (Exception e)
		{
			StsException.outputException("DBInputStream::readObject(Class, Object) failed for object of class " + dbClass,
				e, StsException.WARNING);
			return null;
		}
	}

	public Object readObjectFully(StsDBTypeObject dbClass)
	{
		try
		{
			return dbClass.readObjectFully(this);
		}
		catch (Exception e)
		{
			StsException.outputException("DBInputStream::readObject(Class, Object) failed for object of class " + dbClass,
				e, StsException.WARNING);
			return null;
		}
	}

	public Object readObjectFully(StsDBTypeObject dbClass, Object object)
	{
		try
		{
			return dbClass.readObjectFully(this, object);
		}
		catch (Exception e)
		{
			StsException.outputException("DBInputStream::readObject(Class, Object) failed for object of class " + dbClass,
				e, StsException.WARNING);
			return null;
		}
	}

    public void setPositionEndDB()
    {
        dbFile.setPositionEndDB();
    }
/*
    public boolean moreDataAvailable() throws IOException
    {
        return available() > 0 && !endDBFound;
    }
*/
	public StsDBType getInputDBType(Class c)
	{
		return inputDBClasses.getDBType(c);
	}

	public StsDBType getInputDBType(Object o)
	{
		return inputDBClasses.getDBType(o);
	}

	public StsDBType getInputDBType(int index)
	{
		return inputDBClasses.getDBType(index);
	}

	public StsDBType getInputDBType(String s)
	{
		StsDBType dbClass = inputDBClasses.getDBType(s, null);
		if (dbClass == null)
		{
			throw new RuntimeException("DBInputStream::getInputDBClass(String) - no dbClass found for " + s);
		}
		return dbClass;
	}

	public StsDBObjectTypeList getInputTypeList()
	{
		return inputDBClasses;
	}

	public boolean readBoolean() throws java.io.IOException
	{
		boolean value = super.readBoolean();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public byte readByte() throws java.io.IOException
	{
		byte value = super.readByte();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public char readChar() throws java.io.IOException
	{
		char value = super.readChar();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public double readDouble() throws java.io.IOException
	{
		double value = super.readDouble();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public float readFloat() throws java.io.IOException
	{
		float value = super.readFloat();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public int readInt() throws java.io.IOException
	{
		int value = super.readInt();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public long readLong() throws java.io.IOException
	{
		long value = super.readLong();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public short readShort() throws java.io.IOException
	{
		short value = super.readShort();
		if (debug) System.out.println("Input----->" + value);
		return value;
	}

	public String readUTF() throws java.io.IOException
	{
        try
        {
		    String value = super.readUTF();
            if (debug) System.out.println("Input----->" + value);
            return value;
        }
        catch(Exception e)
        {
            StsException.systemError(this, "readUTF");
            return "bad";
        }
	}

	private class InputDBClasses extends StsDBObjectTypeList
	{
		public InputDBClasses()
		{
			super("inputDBClasses");
		}

		public StsDBTypeClass setDBClass(StsDBTypeClass currentDBClass, StsDBTypeClass dbClass)
		{
			dbClass.setFieldsFromCurrentDBClass(currentDBClass);
			insert(dbClass, dbClass.getIndex());
			return dbClass;
		}
	}
}
