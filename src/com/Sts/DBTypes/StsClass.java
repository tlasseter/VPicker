
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.DB.DBCommand.*;
import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsObjectListPanel;
import com.Sts.UI.StsOkDialog;
import com.Sts.Utilities.*;

import java.awt.event.*;
import java.beans.*;
import java.lang.reflect.*;
import java.util.*;

public class StsClass extends StsSerialize implements StsSerializable, StsTreeObjectI
{
    protected String name;
    protected String userName;
    protected boolean isVisible = true;
    protected StsObject currentObject = null;
    protected int initialSize = 10;
    protected int initialSizeInc = 10;
    protected int lastIndex = -1;

	transient protected StsObjectList list = null;
    transient protected Class c = null;
    transient protected Constructor constructor = null;

    transient public StsFieldBean[] displayFields = null;
    transient public StsFieldBean[] propertyFields = null;
	transient public StsFieldBean[] defaultFields = null;
	transient protected StsObjectPanel objectPanel = null;
	transient public ArrayList<StsClass> subClasses = null;

    transient protected ArrayList actionListeners = new ArrayList();
    transient protected ArrayList propertyChangeListeners = new ArrayList();

    static final boolean debug = false;

   	public StsClass()
    {
    }

    public void initializeListSize(int size, int inc)
    {
        initialSize = size;
        initialSizeInc = inc;
    }

    public void setName(String className)
    {
        this.name = className;
    }

    public void initializeList()
    {
        try
        {
            c = Class.forName(name);
            if(c == null)
            {
                StsException.systemError("StsClass.initializeList() failed. Cannot find class for name: " + name);
                return;
            }
            list = new StsObjectList(initialSize, initialSizeInc, c);
        }
        catch(Exception e)
        {
            StsException.systemError("StsClass.initializeList() failed for " + name);
        }

        try
        {
            constructor = c.getDeclaredConstructor(new Class[0]);
        }
        catch(Exception e)
        {
        	StsException.systemError(this, "initializeList", "Default constructor not found for class " + c);
        }
    }

    public Class getObjectClass() { return c; }

    /** Override as needed in StsClass subclasses to handle initializations when stsClass is constructed */
	public void initializeFields()
	{
	}

    /** For a new project, perform any operations after class has been created. Override in subclasses. */
    public boolean projectInitialize(StsModel model)
    {
		return true;
    }

	/** For a db which has been loaded, perform any operations after class has been created. Override in subclasses. */
	public boolean dbInitialize()
	{
		return true;
    }

    /** db has been loaded, windows are all up and running: final initialization for toolbars, etc. */
    public void finalInitialize()
    {
    }

    public void projectRotationAngleChanged()
    {
        System.out.println("reinitializeRotated() needs to be implemented in class: " + StsToolkit.getSimpleClassname(this));
    }

	/** override in subclasses if display fields name is different than displayFields */
    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null) initializeDisplayFields();
        return displayFields;
    }
	/** override in subclasses if default fields name is different than defaultFields */
	public StsFieldBean[] getDefaultFields()
    {
        if(defaultFields == null) initializeDefaultFields();
        return defaultFields;
    }
	/** override in subclasses if properties fields name is different than propertyFields */
    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields == null)
            initializePropertyFields();
        return propertyFields;
    }

    /** override in subclass if needed */
    public void initializeDisplayFields() { displayFields = null; }
    public void initializeDefaultFields() { defaultFields = null; }
    public void initializePropertyFields() { propertyFields = null; }

    public Object[] getChildren() { return list == null ? null : (Object[])list.copyArrayList(); }

    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
//        else objectPanel.setViewObject(this);
        return objectPanel;
    }

	public void deleteObjectPanel() { objectPanel = null; }
    public void treeObjectSelected() { }

	public int getSize() { return list.getSize(); }

	public StsObject getElement(int index)
    {
        StsObject obj = (StsObject)list.getElement(index);
//        if(obj == null)
//            StsException.systemError("Failed to find instance " + index +
//                " for class: " + c.getName());
        return obj;
    }

    /** returns the object with this index */
    public StsObject getElementWithIndex(int index)
    {
        int size = list.getSize();
        if(size == 0) return null;
        StsObject obj = (StsObject)list.getElement(index);
        if(obj != null && obj.getIndex() == index) return obj;
		int nStart = size-1;
        //StsObject lastObject = (StsObject)list.getElement(nStart);
        //if(lastObject.getIndex() < index)
        //    return null;
        //        int nStart = Math.min(index, size-1);
        for(int n = nStart; n >= 0; n--)
        {
            obj = (StsObject)list.getElement(n);
            if(obj.getIndex() == index)
            {
                return obj;
            }
        }
        return null;
    }
    public String getName() { return name; }
    public void setIsVisible(boolean isVisible)
    {
        this.isVisible = isVisible;
    }

    public boolean getIsVisible() { return isVisible; }

	public StsObject getFirst() { return (StsObject)list.getFirst(); }
    public StsObject getNext() { return (StsObject)list.getNext(); }
	public StsObject getLast() { return (StsObject)list.getLast(); }
    public Enumeration elements() { return list.elements(); }
	public boolean contains(StsObject object) { return list.contains(object); }

    public StsObject getCurrentObject() { return currentObject; }

    public boolean setCurrentObject(StsObject object)
    {
        if(currentObject == object) return false;
 //       currentModel.viewObjectChangedAndRepaint(object);
        currentObject = object;
//        setViewObject(object);
        return true;
    }

    /** Override in subclasses which wish to inform views that currentObject has changed. */
/*
    public void setViewObject(StsObject object)
    {
    }
*/
    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if(!isVisible) return;

        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsObject object = (StsObject)iter.next();
            object.display(glPanel3d);
        }
    }

    public void add(StsObject obj) throws StsException
    {
        try
        {
            if(obj == null) throw new StsException(StsException.WARNING,
                                "StsClass.add(StsObject)", "Tried to add null object");

			// It is possible that we have a reference to an object before the object itself has been instantiated.
			// So insert it in the list based on its index number.
			if(obj.getIndex() != -1)
			{
				if(checkInsertInList(obj)) return;
				StsException.systemError("Object of stsClass " + getName() + " already in instance list; developer you are adding it twice!");
				return;
			}

            list.add(obj);
            obj.setIndex(++lastIndex);
            StsInstanceAddCmd cmd = new StsInstanceAddCmd(obj);
            obj.getCurrentModel().addTransactionCmd("instanceListAdd", cmd);
            setCurrentObject(obj);
            fireActionEvent(obj);
        }
        catch(Exception e)
        {
            StsException.outputException("StsClass.add() failed.", e, StsException.WARNING);
        }
    }

    public void addActionListener(ActionListener listener)
    {
        actionListeners.add(listener);
    }

    public void removeActionListener(ActionListener listener)
    {
        actionListeners.remove(listener);
    }

    public void fireActionEvent(Object obj)
    {
        int nActionListeners = actionListeners.size();
        if(nActionListeners == 0) return;

        ActionEvent actionEvent = new ActionEvent(obj, 0, "add");
        for (int n = 0; n < nActionListeners; n++)
        {
            ActionListener listener = (ActionListener) actionListeners.get(n);
            listener.actionPerformed(actionEvent);
        }
    }

    public void addPropertyChangeListener(PropertyChangeListener listener)
    {
        propertyChangeListeners.add(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener)
    {
        propertyChangeListeners.remove(listener);
    }

    public void firePropertyChanged(PropertyChangeEvent changeEvent)
    {
        int nChangeListeners = propertyChangeListeners.size();
        if(nChangeListeners == 0) return;

        for (int n = 0; n < nChangeListeners; n++)
        {
            PropertyChangeListener listener = (PropertyChangeListener)propertyChangeListeners.get(n);
            listener.propertyChange(changeEvent);
        }
    }

    public void add(StsObject[] objs) throws StsException
    {
        if (objs==null)
        	throw new StsException(StsException.FATAL, "StsClass.add(StsObject[])",
                    "Tried to add null object array");
        for (int i=0; i<objs.length; i++) add(objs[i]);
    }

    public StsDBCommand createInstanceChgCmd(StsObject obj)
    {
        return new StsInstanceChgCmd(obj);
    }

    public void delete(StsObject object)
    {
        if(object == null) return;
        int deleteIndex = -1;
        int size = getSize();
        if (size < 1) return;
        int index = object.getIndex();
        if(getElement(index) == object)
            deleteIndex = index;
        else
        {
            int indexStart = size-1;
    //        int indexStart = Math.min(index, size-1);
            for(int n = indexStart; n >= 0; n--)
            {
                StsObject matchObject = getElement(n);
                if(matchObject == object)
                {
                    deleteIndex = n;

                }
            }
        }
        if(deleteIndex == -1) return;
        delete(object, deleteIndex);
        currentObject = currentModel.deleteCursor3dObject(object);
	    setCurrentObject(currentObject);
    }

    public boolean deleteObjectWithIndex(int index)
    {
        StsObject object;
        int indexStart = 0;
        int size = getSize();
        if (size < 1) return false;
        if(index < size)
        {
            object = getElement(index);
            if(object.getIndex() == index)
            {
                delete(object, index);
                return true;
            }
            indexStart = index;
        }
        else
            indexStart = size-1;

        for(int n = indexStart; n >= 0; n--)
        {
            StsObject matchObject = getElement(n);
            if(matchObject.getIndex() == index)
            {
                delete(matchObject, n);
                return true;
            }
        }
        return false;
    }

	public void rundelete(StsObject s, ActionListener a)
	  {
		  final ActionListener aa = a;
		  final StsObject ss = s;
		  Runnable runnable = new Runnable()
		  {
			  public void run()
			  {
				  aa.actionPerformed(new ActionEvent(ss, 0, "delete"));
			  }
		  };
		  StsToolkit.runWaitOnEventThread(runnable);
	  }


    private void delete(StsObject object, int sequentialIndex)
    {
        list.delete(sequentialIndex);
		StsInstanceDeleteCmd cmd = constructInstanceDeleteCmd(object);
        object.getCurrentModel().addTransactionCmd("deleteInstanceObj", cmd);
        object.setIndex(-1);

        for(int n = 0; n < actionListeners.size(); n++)
        {
            ActionListener listener = (ActionListener)actionListeners.get(n);
			rundelete(object,listener);
			// jbw thread bind
            //listener.actionPerformed(new ActionEvent(object, 0, "delete"));
        }
        if(list.getSize() == 0)
        {
            setCurrentObject(null);
            currentModel.deleteCursor3dTextures(object);
        }
        else
            setCurrentObject((StsObject)list.getLast());

        currentModel.refreshObjectPanel(object);
    }

    public boolean hasInterface(Class classInterface)
    {
        Class[] interfaces = getClass().getInterfaces();
        if(interfaces == null) return false;
        for(int i = 0; i < interfaces.length; i++)
            if(interfaces[i] == classInterface) return true;
        return false;
    }

    public void selected(StsObject object)
    {
        for(int n = 0; n < actionListeners.size(); n++)
        {
            ActionListener listener = (ActionListener)actionListeners.get(n);
            listener.actionPerformed(new ActionEvent(object, 0, "selected"));
        }
    }

    public StsInstanceDeleteCmd constructInstanceDeleteCmd(StsObject obj)
    {
        int index = obj.getIndex();
        if (index < 0)
        {
            StsException.systemError("StsClass.InstanceDeleteCmd.constructor() failed. Index is " + index + " for object " + obj.getName());
            return null;
        }
        return new StsInstanceDeleteCmd(obj);
    }

    public StsObject getEmptyStsClassObject(int index)
    {
		StsObject object = null;
        try
        {
            object = (StsObject)constructor.newInstance(new Object[0]);
            object.setIndex(index);
			// checkListSequence(object); // this is a hack:  try to fix!
            return object;
        }
        catch(Exception e)
        {
            StsException.outputException("StsClass.getEmptyStsClassObject() failed for class " + name,
                e, StsException.WARNING);
            return null;
        }
    }

    private boolean checkInsertInList(StsObject object)
    {
		int index = object.getIndex();
		int nElements = getSize();
		if(nElements == 0 || nElements == index)
		{
			list.add(object);
			lastIndex = index;
			return true;
		}
		else
		{
			int nextIndex = -1;
			for(int n = 0; n < nElements; n++)
			{
				int prevIndex =  nextIndex;
				nextIndex = ((StsObject)list.getElement(n)).getIndex();
				if(index > prevIndex)
				{
					if(index == nextIndex)
					{
						StsException.systemError("System error: duplicate indexes for list in StsClass " + getName());
						return false;
					}
					list.insertBefore(n, object);
					return true;
				}
			}
			list.add(object);
			lastIndex = index;
			return true;
		}
	}

    private void checkListSequence(StsObject object)
    {
        StsObject objectAtIndex;
		int index = object.getIndex();
        int nElements = list.getSize();
        if(index < nElements)
        {
            objectAtIndex = (StsObject)list.getElement(index);
            if(object == objectAtIndex) return;
            int indexAt = objectAtIndex.getIndex();
            if(index < indexAt)
                list.insertBefore(indexAt, object);
            else
                list.insertAfter(indexAt, object);
        }
        list.insertAt(index, object);
	}

    public void preAllocate(int numInstances)
    {
    	try
        {
	    	for( int j=0; j<numInstances; j++ )
    	    	constructor.newInstance(new Object[0]);
        }
        catch(Exception e)
        {
        	System.out.println("Exception in preAllocate()\n" + e);
        }
    }

    /** can this calss export all data in the view **/
    public boolean canExportView() { return false; }
    public int exportView()
    {
        StsMessageFiles.infoMessage("Objects of " + name + " cannot export objects in view.");
        return 0;
    }

    /** return an array of names */
    public String[] getNames()
    {
        if (list==null) return new String[0];
        return list.getNames();
    }

    public StsList getList() { return list; }
    public StsObjectList getStsObjectList() { return list; }

    /** get a list with objects in the order of name array */
    public StsList getObjectList(String[] names)
    {
        if (list==null || names == null) return null;
        int nItems = getSize();
        if (nItems==0) return null;
        int nObjects = Math.min(nItems, names.length);
        StsList list = new StsList(nObjects);
        for (int i=0; i<nObjects; i++)
        {
            StsObject obj = getElement(names[i]);
            if (obj!=null) list.add(obj);
        }
        return list;
    }

	public Object[] getTrimmedList() { return list.getTrimmedList(); }

    /** get a list with objects excluding these objects */
    public StsObjectList getObjectListExcluding(StsObjectList excludedObjects)
    {
        if (list==null || excludedObjects == null) return null;
        int nItems = getSize();
        if (nItems==0) return null;
        StsObjectList list = new StsObjectList(nItems);
        for (int i=0; i<nItems; i++)
        {
            StsObject obj = getElement(i);
            StsObject excludedObject = excludedObjects.getObjectWithName(obj.getName());
            if(excludedObject == null) list.add(obj);
        }
        return list;
    }
    /** get a list with objects of this type excluding these objects */
    public StsObjectList getObjectListOfTypeExcluding(byte type, StsObjectList excludedObjects)
    {
        if (list == null) return null;

        int nItems = getSize();
        if (nItems==0) return null;
        StsObjectList list = new StsObjectList(nItems);
        StsObject excludedObject = null;
        for (int i=0; i<nItems; i++)
        {
            StsMainObject obj = (StsMainObject)getElement(i);
            if(obj.getType() != type) continue;
            if(excludedObjects == null)
                list.add(obj);
            else
            {
                excludedObject = excludedObjects.getObjectWithName(obj.getName());
                if (excludedObject == null) list.add(obj);
            }
        }
        return list;
    }

    /** return an object from its name */
    public StsObject getElement(String name)
    {
        if (list==null) return null;
        return list.getElement(name);
    }

    /** return an array of objects */
    public StsObject[] getElements() { return list.getElements(); }

	public Object getCastObjectList()
	{
        int nElements = getSize();
        if(nElements == 0)
            return Array.newInstance(c, 0);
		else
            return StsMath.arraycopy(list.getList(), nElements, c);
	}

    public Object getCastObjectList(Class c)
    {
        int nElements = getSize();
        if(nElements == 0)
            return Array.newInstance(c, 0);
        else
            return StsMath.arraycopy(list.getList(), nElements, c);
    }

	/** Reorder the current objectInstanceList by these names.
	 *  Objects not in this name list will be put at the end.
	 */
	public boolean reorderListByNames(String[] names)
	{
		StsObjectList newList = list.reorderListByNames(names);
		if(newList == null) return false;
		list = newList;
		return true;
	}


	/** Reorder the current objectInstanceList like these objects.
	 *  Objects not in this name list will be put at the end.
	 */
	public boolean reorderList(StsObject[] newObjects)
	{
		StsObjectList newList = list.reorderList(newObjects);
		if(newList == null) return false;
		list = newList;
		return true;
	}

	/** If a single argument, make an array of length 1 for compatibility */
 	public void forEach(String methodName, Object arg)
    {
        if(list == null) return;
    	list.forEach(methodName, new Object[]{arg});
    }

    /** For each element in the list, apply this method/args. */
    public void forEach(String methodName, Object[] args)
    {
        if(list == null) return;
        list.forEach(methodName, args);
    }

    public void forEach(String methodName)
    {
        if(list == null) return;
        list.forEach(methodName, null);
    }

	public void forEachInvoke(Method method, Object[] args)
	{
        if(list == null) return;
		list.forEachInvoke(method, args);
	}

    public StsObjectList.ObjectIterator getObjectIterator()
    {
        return list.getObjectIterator(false);
    }

    public StsObjectList.ObjectIterator getObjectIterator(boolean visibleOnly)
    {
        return list.getObjectIterator(visibleOnly);
    }

    public StsObjectList.ObjectIterator getObjectIterator(boolean visibleOnly, byte type)
    {
        return list.getObjectIterator(visibleOnly, type);
    }

    public StsObjectList.ObjectIterator getObjectIterator(byte type)
    {
        return list.getObjectIterator(type);
    }

    public StsObjectList.ObjectIterator getVisibleObjectIterator()
    {
        return list.getVisibleObjectIterator();
    }

    public StsObjectList.ObjectIterator getObjectOfTypeIterator(byte type)
    {
        return list.getObjectOfTypeIterator(type);
    }

    public StsObjectList.ObjectIterator getVisibleObjectOfTypeIterator(byte type)
    {
        return list.getVisibleObjectOfTypeIterator(type);
    }

    public StsObject[] getObjectList()
    {
        return list.getObjectList();
    }

    public StsMainObject[] getVisibleObjectList()
    {
        return list.getVisibleObjectList();
    }

    public StsMainObject[] getObjectListOfType(byte type)
    {
        return list.getObjectListOfType(type);
    }

    public Object[] getObjectListOrderedByName()
    {
        return list.getObjectListOrderedByName();
    }

    public StsMainObject[] getVisibleObjectListOfType(byte type)
    {
        return list.getVisibleObjectListOfType(type);
    }

    public StsObject getObjectWithName(String name)
    {
        return list.getObjectWithName(name);
    }

    public StsMainObject getObjectOfType(byte type, int index)
    {
        return list.getObjectOfType(type, index);
    }

    public boolean hasObjectsOfType(byte type)
    {
        return list.hasObjectsOfType(type);
    }

    public boolean hasObjects()
    {
        return list.getSize() > 0;
    }

    public int getNObjectsOfType(byte type)
    {
        return list.getNObjectsOfType(type);
    }
/*
    public void delete(StsClass stsObject)
    {
        for(int n = 0; n < actionListeners.size(); n++)
        {
            ActionListener listener = (ActionListener)actionListeners.get(n);
            listener.actionPerformed(new ActionEvent(stsObject, 0, "delete"));
        }
    }
*/
    public void deleteAll()
    {
        // delete each entry
        int nObjects = getSize();
        for(int n = nObjects-1; n >= 0 ; n--)
        {
            StsObject object = getElement(n);
            object.delete();
        }
        // remove all entries from list
        this.list.deleteAll();
    }

   public void deleteNonPersistentObjects()
    {
        // delete each entry
        int nObjects = getSize();
        for(int n = nObjects-1; n >= 0 ; n--)
        {
            StsObject object = getElement(n);
            if(object.isPersistent()) continue;
            object.delete();
        }
    }

    public StsObject setPreviousObject()
    {
        StsObject prevObject = getPreviousObject(currentObject);
        if(prevObject != null) setCurrentObject(prevObject);
        return prevObject;
    }

    public StsObject setNextObject()
    {
        StsObject nextObject = getNextObject(currentObject);
        if(nextObject != null) setCurrentObject(nextObject);
        return nextObject;
    }

    /** Allow for the fact that an object may have been deleted so that the index
     * of this object is greater than its current sequential position. So search
     * down from the index position until you have a match and return previous one.
     */
    public StsObject getPreviousObject(StsObject object)
    {
        if(object == null) return null;
        int size = getSize();
        if (size <= 1) return null;
 //       int index = object.index();
 //       if(index == 0) return null;

        for(int n = size-1; n >= 0; n--)
        {
            if(getElement(n) == object)
			{
				if (n > 0)
					return getElement(n - 1);
				else
					return getElement(size-1);
			}
        }
        return null;
    }

    /** Allow for the fact that an object may have been deleted so that the index
     * of this object is greater than its current sequential position. So search
     * down from the index position until you have a match and return the next one.
     */
    public StsObject getNextObject(StsObject object)
    {
		if(object == null) return null;
	    int size = getSize();
	    if (size <= 1) return null;
//       int index = object.index();
//       if(index == 0) return null;

	    for(int n = 0; n < size; n++)
	    {
			if (getElement(n) == object)
			{
				if (n < size - 1)
					return getElement(n + 1);
				else
					return getElement(0);
			}
		}
		return null;
    }

    // need signatures for each win3d type because reflection in StsMethod
    // cant find a general signature

    public void toggleOn(StsWin3dBase win3d)
    {
        if(debug) System.out.println("toggleOn called.");
        win3d.getCursor3d().toggleOn(getCurrentObject());
        win3d.repaint();
    }

    public void toggleOff(StsWin3dBase win3d)
    {
        if(debug) System.out.println("toggleOff called.");
        win3d.getCursor3d().toggleOff(getCurrentObject());
        win3d.repaint();
    }

    public void toggleOn(StsWin3d win3d)
    {
        if(debug) System.out.println("toggleOn called.");
        win3d.getCursor3d().toggleOn(getCurrentObject());
        win3d.repaint();
    }

    public void toggleOff(StsWin3d win3d)
    {
        if(debug) System.out.println("toggleOff called.");
        win3d.getCursor3d().toggleOff(getCurrentObject());
        win3d.repaint();
    }

    public void toggleOn(StsWin3dFull win3d)
    {
        if(debug) System.out.println("toggleOn called.");
        win3d.getCursor3d().toggleOn(currentObject);
        win3d.repaint();
    }

    public void toggleOff(StsWin3dFull win3d)
    {
        if(debug) System.out.println("toggleOff called.");
        win3d.getCursor3d().toggleOff(getCurrentObject());
        win3d.repaint();
    }

    public void close() { }

    public boolean anyDependencies() { return true; }
    public boolean export() { return false; }
    public boolean launch() { return false; }
    public boolean canExport() { return false; }
    public boolean canLaunch() { return false; }
    public void popupPropertyPanel() { return; }

    /** implement in any subclasses which use displayLists */

    public void deleteDisplayLists()
    {
    }

    public StsFieldBean getBeanWithName(String beanName)
    {
        if(propertyFields != null)
        {
            StsFieldBean bean = StsFieldBean.getBeanWithFieldName(propertyFields, beanName);
            if(bean != null) return bean;
        }
        if(displayFields != null)
        {
            StsFieldBean bean = StsFieldBean.getBeanWithFieldName(displayFields, beanName);
            if(bean != null) return bean;
        }
		if(defaultFields != null)
		{
			StsFieldBean bean = StsFieldBean.getBeanWithFieldName(defaultFields, beanName);
			if(bean != null) return bean;
        }
        return null;
    }

    public Class[] getViewClasses()
    {
        if(getSize() == 0) return new Class[0];
        if(!StsViewable.class.isAssignableFrom(c)) return new Class[0];
        if(currentObject == null) return new Class[0];
        StsViewable viewable = (StsViewable)currentObject;
        return viewable.getViewClasses();
    }

	public ArrayList<StsClass> getSubClasses()
	{
		if(subClasses != null) return subClasses;
		return null;
	}

	public void addSubClass(StsClass subClassStsClass)
	{
		if(subClasses == null) subClasses = new ArrayList<StsClass>();
		subClasses.add(subClassStsClass);
	}

    public String toString() { return name; }


    // Can only be used in viewer since the beans are static. The object panel is not shown when in viewer only mode so not a problem.
    public void popupClassListPanel()
    {
		try
        {
            StsObjectListPanel box = new StsObjectListPanel(currentModel, this);
            if(box.initialize())
            {
                box.gbc.fill = box.gbc.HORIZONTAL;
                new StsOkDialog(currentModel.win3d, box, getName() + " List", false);
            }
            else
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING, "No Objects of this type in project.", true);
            }
        }
        catch(Exception ex)
        {
            StsException.outputException("Failed to construct class object list panel for " + this.getName(), ex, StsException.WARNING);
        }
    }

    public String getUserName()
    {
        if(userName == null)
            return name;
        else
            return userName;
    }

/*
	static public void initColors(StsFieldBean[] fields)
	{
		StsSpectrum spectrum = currentModel.getSpectrum("Basic");
		StsColorListFieldBean colorList = (StsColorListFieldBean) StsFieldBean.getBeanWithFieldName(fields, "color");
		if (colorList != null) colorList.setListItems(spectrum);
	}*/
}
