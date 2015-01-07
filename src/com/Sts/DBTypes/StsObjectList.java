
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsObjectList extends StsList implements StsSerializable
{
 //   protected String className = null;
 //   transient Class c = null;

    public StsObjectList()
    {
    }

    public StsObjectList(int initSize, int incSize)
    {
        super(initSize, incSize);
    }

    public StsObjectList(int initSize, int incSize, Class c)
    {
        super(initSize, incSize);
    //    className = c.getName();
    //    this.c = c;
    }

   	public StsObjectList(int initSize)
    {
        this(initSize, 1);
    }

	public StsObjectList(StsObject[] stsObjects, int incSize)
	{
		super(stsObjects.length, incSize);
		add(stsObjects);
    }

    private StsObjectList(StsList list)
    {
        super(list);
    }

    static public StsObjectList constructor(StsList list)
    {
        if(list instanceof StsObjectList)
            return (StsObjectList)list;
        
        // verify that all list objects are subclassed from StsObject
        for(Object object : list.getList())
        {
            if(!(object instanceof StsObject))
            {
                StsException.systemError(StsObjectList.class, "constructor", "Object " + object.toString() + " is not a subclass instance of StsObject");
                return null;
            }
        }
        return new StsObjectList(list);
    }

    public boolean initialize(StsModel model)
    {
        return true;
//        return setClassFromName();
    }
/*
    private boolean setClassFromName()
    {
        try
        {
            if(className != null) c = Class.forName(className);
            return true;
        }
        catch(Exception e)
        {
            StsException.systemError("StsObjectList.classInitialize() failed." +
                " Could not find class: " + className);
            return false;
        }
    }
*/
    public int add(StsObject obj)
    {
    /*
        Class clazz = obj.getClass();

        if(className == null)
        {
            c = obj.getClass();
            className = c.getName();
        }
        else if(clazz!= c)
        {
            StsException.systemError("StsObjectList.add() failed." +
                " Object: " + obj.getName() + " is not an instance of class " + className);
            return StsList.NO_MATCH;
        }
     */
    // DEBUG!! remove when finished
//        if(obj instanceof StsWellPlanSet)
//            System.out.println("Adding StsWellPlanSet with index " + obj.index());
        return super.add(obj);
    }

    public void add(StsObject[] objs) { super.add(objs); }
/*
    public int delete(StsObject obj)
    {
        return super.delete((Object)obj);
    }
*/
 	public void addList(StsObjectList extraList)
    {
        super.addList((StsList)extraList);
    }

	public int insertBefore(int index, StsObject newObj)
    {
        return super.insertBefore(index, (Object)newObj);
    }

	public int insertBefore(StsObject obj, StsObject newObj)
    {
        return super.insertBefore((Object)obj, (Object)newObj);
    }

    public int insertAfter(int index, StsObject newObj)
    {
        return super.insertAfter(index, (Object)newObj);
    }

    public void unInsert(StsObject obj)
    {
        super.unInsert((Object)obj);
    }

	public int insertAfter(StsObject obj, StsObject newObj)
    {
        return super.insertAfter((Object)obj, (Object)newObj);
    }

    public boolean contains(StsObject object)
    {
        return super.contains((Object)object);
    }

    /** return an array of names */
    public String[] getNames()
    {
        int nItems = getSize();
        if (nItems==0) return null;
        String[] names = new String[nItems];
        for (int i=0; i<nItems; i++)
        {
            StsObject obj = (StsObject)list[i];
            names[i] = obj.getName();
        }
        return names;
    }

    /** return an object from its name */
    public StsObject getElement(String name)
    {
        if (list==null || name==null) return null;
        int nItems = getSize();
        if (nItems==0) return null;
        for (int i=0; i<nItems; i++)
        {
            StsObject obj = (StsObject)list[i];
            if (obj!=null) if (name.equals(obj.getName())) return obj;
        }
        return null;  // no match
    }

    /** return an array of objects */
    public StsObject[] getElements()
    {
        if (list == null) return new StsObject[0];
        int nItems = getSize();
        if (nItems == 0) return new StsObject[0];
        StsObject[] elements = new StsObject[nItems];
        for (int i=0; i<nItems; i++) elements[i] = (StsObject)getElement(i);
        return elements;
    }

	/** get a copy which can be edited without destroying original list */
	public StsObject[] getListCopy() { return getElements(); }

	/** get a cast list copy which can be edited without destroying original list */
	public Object getCastListCopy()
	{
		return StsMath.arraycopy(list, size);
	}

    /** temp method until sorting problem fixed */
    public void checkSortIndexes()
    {
        if (list == null)return;
        int nElements = getSize();
        if (nElements <= 1)return;
        boolean sorted = false;
        while (!sorted)
        {
            sorted = true;
            StsObject nextObject = (StsObject)list[0];
            for (int n = 1; n < nElements; n++)
            {
                StsObject prevObject = nextObject;
                nextObject = (StsObject)list[n];
                if (prevObject.getIndex() == -1)
                {
                    StsException.systemError("StsObjectList.checkSortIndexes() error: index = -1 for object " +
                                             prevObject.getName());
                    continue;
                }
                if (nextObject.getIndex() == -1)
                {
                    StsException.systemError("StsObjectList.checkSortIndexes() error: index = -1 for object " +
                                             nextObject.getName());
                    continue;
                }
                if (prevObject.getIndex() > nextObject.getIndex())
                {
                    list[n-1] = nextObject;
                    list[n] = prevObject;
                    sorted = false;
                }
            }
        }
    }

    public ObjectIterator getObjectIterator()
    {
        return new ObjectIterator(false);
    }

    public ObjectIterator getObjectIterator(boolean visibleOnly)
    {
        return new ObjectIterator(visibleOnly);
    }

    public ObjectIterator getObjectIterator(boolean visibleOnly, byte type)
    {
        return new ObjectIterator(visibleOnly, type);
    }

    public ObjectIterator getObjectIterator(byte type)
    {
        return new ObjectIterator(false, type);
    }

    public ObjectIterator getVisibleObjectIterator()
    {
        return new ObjectIterator(true);
    }

    public ObjectIterator getObjectOfTypeIterator(byte type)
    {
        return new ObjectIterator(false, type);
    }

    public ObjectIterator getVisibleObjectOfTypeIterator(byte type)
    {
        return new ObjectIterator(true, type);
    }

    public StsObject[] getObjectList()
    {
        return getElements();
    }

    public StsMainObject[] getVisibleObjectList()
    {
        StsMainObject[] objectsList = new StsMainObject[size];
        int nn = 0;
        for(int n = 0; n < size; n++)
        {
            StsMainObject object = (StsMainObject)list[n];
            if(object.getIsVisible()) objectsList[nn++] = object;
        }
        return (StsMainObject[])StsMath.arraycopy(objectsList, nn, StsMainObject.class);
    }

    public StsMainObject[] getObjectListOfType(byte type)
    {
        StsMainObject[] objectsList = new StsMainObject[size];
        int nn = 0;
        for(int n = 0; n < size; n++)
        {
            StsMainObject object = (StsMainObject)list[n];
            if(object.isType(type)) objectsList[nn++] = object;
        }
        return (StsMainObject[])StsMath.arraycopy(objectsList, nn, StsMainObject.class);
    }

    public Object[] getObjectListOrderedByName()
    {
        ArrayList objectList = new ArrayList(size);
        for(int n = 0; n < size; n++)
            objectList.add(list[n]);

        Comparator comparator = StsMainObject.getNameComparator();
        Collections.sort(objectList, comparator);
        return objectList.toArray();
    }

    public StsMainObject[] getVisibleObjectListOfType(byte type)
    {
        StsMainObject[] objectsList = new StsMainObject[size];
        int nn = 0;
        for(int n = 0; n < size; n++)
        {
            StsMainObject object = (StsMainObject)list[n];
            if(object.getIsVisible() && object.isType(type)) objectsList[nn++] = object;
        }
        return (StsMainObject[])StsMath.arraycopy(objectsList, nn, StsMainObject.class);
    }

    public StsObject getObjectWithName(String name)
    {
        for(int n = 0; n < size; n++)
        {
            StsObject object = (StsObject)list[n];
            if(object == null) continue;
            String objectName = object.getName();
            if(objectName == null) continue;
            if(objectName.equals(name)) return object;
        }
        return null;
    }

    public StsMainObject getObjectOfType(byte type, int index)
    {
        int nn = 0;
        for(int n = 0; n < size; n++)
        {
            StsMainObject object = (StsMainObject)list[n];
            if(object.isType(type) && nn++ == index) return object;
        }
        return null;
    }

    public boolean hasObjectsOfType(byte type)
    {
        ObjectIterator iterator = new ObjectIterator();
        while(iterator.hasNext())
        {
            StsMainObject object = (StsMainObject)iterator.next();
            if(object.isType(type)) return true;
        }
        return false;
    }

    public boolean hasObjects()
    {
        return size > 0;
    }

    public int getNObjectsOfType(byte type)
    {
        ObjectIterator iterator = new ObjectIterator();
        int nObjects = 0;
        while(iterator.hasNext())
        {
            StsMainObject object = (StsMainObject)iterator.next();
            if(object.isType(type)) nObjects++;
        }
        return nObjects;
    }

    public class ObjectIterator implements Iterator
    {
        int nObjects = 0;
        int n = 0;
        StsMainObject next = null;
        boolean visibleOnly = false;
        byte type = -128;

        ObjectIterator()
        {
            nObjects = size;
            if(nObjects > 0) setNext();
        }

        ObjectIterator(boolean visibleOnly)
        {

            nObjects = size;
            this.visibleOnly = visibleOnly;
            if(nObjects > 0) setNext();
        }

        ObjectIterator(boolean visibleOnly, byte type)
        {
            nObjects = size;
            this.visibleOnly = visibleOnly;
            if(type < 0 || type == -128)
                StsException.systemError("ObjectIterator constructor()." +
                    " Illegal type: -128. Must be from -127 to 127.");
            this.type = type;
            if(nObjects > 0) setNext();
        }

        private void setNext()
        {
            if(visibleOnly)
            {
                if(type == -128)
                {
                    while(n < size)
                    {
                        next = (StsMainObject)list[n++];
                       if(next.getIsVisible()) return;
                    }
                }
                else
                {
                    while(n < size)
                    {
                        next = (StsMainObject)list[n++];
                        if(next.getIsVisible() && next.isType(type)) return;
                    }
                }
            }
            else
            {
                if(type != -128)
                {
                    while(n < size)
                    {
                        next = (StsMainObject)list[n++];
                        if(next.isType(type)) return;
                    }
                }
                else
                {
                    while(n < size)
                    {
                        next = (StsMainObject)list[n++];
                        return;
                    }
                }
            }
            next = null;
        }

        public boolean hasNext()
        {
            return next != null;
        }

        public Object next()
        {
            StsMainObject current = next;
            setNext();
            return current;
        }

        public void remove()
        {
        }

        public int getSize()
        {
            StsMainObject object;
            int n = 0, nn = 0;

            if(visibleOnly)
            {
                if(type == -128)
                    while(n < size)
                    {
                        object = (StsMainObject)list[n++];
                        if(object.getIsVisible()) nn++;
                    }
                else
                    while(n < size)
                    {
                        object = (StsMainObject)list[n++];
                        if(object.getIsVisible() && object.isType(type)) nn++;
                    }
            }
            else
            {
                if(type != -128)
                    while(n < size)
                    {
                        object = (StsMainObject)list[n++];
                        if(object.isType(type)) nn++;
                    }
                else
                    return nObjects;
            }
            return nn;
        }
    }

    public StsObjectRefList convertListToRefList(StsModel model, String fieldName, StsObject parentObj)
    {
        for(int n = 0; n < size; n++)
        {
            StsObject object = (StsObject)list[n];
            if(!object.isPersistent()) object.addToModel();
        }
        return StsObjectRefList.convertListToRefList(this, fieldName, parentObj);
    }

	public StsObjectList reorderListByNames(String[] names)
	{
		StsObjectList newList = new StsObjectList(size, incSize);
		int nNames = names.length;
		for(int n = 0; n < nNames; n++)
		{
			StsObject object = getObjectWithName(names[n]);
			if(object == null) continue;
			newList.add(object);
			delete(object);
		}
		int nRemainingObjects = this.getSize();
		for(int n = 0; n < nRemainingObjects; n++)
			newList.add(getElement(n));
		return newList;
	}

	public StsObjectList reorderList(StsObject[] objects)
	{
		StsObjectList newList = new StsObjectList(size, incSize);
		newList.add(objects);
		return newList;
	}

    public StsObjectList deepCopy()
    {
        StsObject[] newObjects = new StsObject[size];
        for(int n = 0; n < size; n++)
        {
            Object newObject = ((StsObject)list[n]).clone();
            newObjects[n] = (StsObject)newObject;
        }
        return new StsObjectList(newObjects, incSize);
    }

/*
	public void write(DBOutputStream out)
	{
   		// write out the size and increment
        try
        {
            out.writeInt(size);
            if(size == 0) return;

            if(className == null) out.writeUTF("none");
            else                  out.writeUTF(className);

            // write out all of the non-null entries
            for( int i=0; i<size; i++ )
            {
                StsObject obj = (StsObject)list[i];
                if(obj != null) out.write(obj, c);
            }
            if(Main.isDbDebug) System.out.println("     write StsObjectList: " + className + "[" + size + "]");
        }
        catch(Exception e)
        {
        	System.out.println("Exception in StsObjectList.write()\n" + e);
        }
    }

    public void read(DBInputStream in)
    {
        int n = -1;
        try
        {

            // read the size and class
            int size = in.readInt();
			className = (String) in.readUTF();
            c = Class.forName(className);

            // read the entries and associate with the array
            for(n = 0; n < size; n++ )
                in.read(list[n], c);
        }
        catch(Exception e)
        {
            String name;
            if(className != null) name = className + "[" + n + "]";
            else             name = "null";
            StsException.outputException("Error in StsClass.read() for instance " + name,
                e, StsException.WARNING);
        }
    }
*/
}
