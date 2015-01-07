
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Beans;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.lang.reflect.*;
import java.util.*;

public class StsComboBoxFieldBean extends StsFieldBean implements PopupMenuListener,ActionListener,ItemListener
{
    JComboBox comboBox = new JComboBox();
    Object[] items;
    Method getItems;
	String listFieldName;

    public StsComboBoxFieldBean()
    {
    }

    public StsComboBoxFieldBean(Object instance, String fieldName, String fieldLabel, Object[] items)
    {
        super();
        initialize(instance, fieldName, fieldLabel, items);
    }

    public StsComboBoxFieldBean(Object instance, String fieldName, String fieldLabel)
    {
        super();
        initialize(instance, fieldName, fieldLabel, (Object[])null);
    }

    public StsComboBoxFieldBean(Object instance, String fieldName, Object[] items)
    {
        super();
        initialize(instance, fieldName, null, items);
    }

    public StsComboBoxFieldBean(Object instance, String fieldName, String fieldLabel, String listFieldName)
    {
        super();
        initialize(instance, fieldName, fieldLabel);
        initializeGetItems(beanObjectClass, listFieldName);
        setItemsFromBeanObject();
    }

    public StsComboBoxFieldBean(Class c, String fieldName, String fieldLabel, Object[] items)
    {
        super();
        classInitialize(c, fieldName, fieldLabel, items);
    }

    public StsComboBoxFieldBean(Class c, String fieldName, String fieldLabel)
    {
        super();
        classInitialize(c, fieldName, fieldLabel);
    }

    public StsComboBoxFieldBean(Class c, String fieldName, String fieldLabel, String listFieldName)
    {
        super();
        classInitialize(c, fieldName, fieldLabel);
        initializeGetItems(c, listFieldName);
    }

    public void classInitialize(Class c, String fieldName, String fieldLabel, Object[] items)
    {
        super.classInitialize(c, fieldName, true, fieldLabel);
        setListItems(items);
//		comboBox.setEnabled(editable);
        comboBox.setEditable(false);
		layoutBean();
        addActionListener();
        setLightWeightPopupEnabled(false);
    }

    public void classInitialize(Class c, String fieldName, String fieldLabel)
    {
        super.classInitialize(c, fieldName, true, fieldLabel);
//		comboBox.setEnabled(editable);
        comboBox.setEditable(false);
		layoutBean();
        addActionListener();
        setLightWeightPopupEnabled(false);
    }

    public void initialize(Object instance, String fieldName, String fieldLabel, Object[] items)
    {
        super.initialize(instance, fieldName, true, fieldLabel);
		setListItems(items);
        comboBox.setEditable(editable);
        layoutBean();
        addActionListener();
        setLightWeightPopupEnabled(false);
    }

    public void initialize(Object instance, String fieldName, String fieldLabel)
    {
        super.initialize(instance, fieldName, true, fieldLabel);
        comboBox.setEditable(editable);
        layoutBean();
        addActionListener();
        setLightWeightPopupEnabled(false);
    }

    public void initialize(Object instance, String fieldName, String fieldLabel, String listFieldName)
    {
        super.initialize(instance, fieldName, true, fieldLabel);
        initializeGetItems(beanObjectClass, listFieldName);
        comboBox.setEditable(editable);
        layoutBean();
        addActionListener();
        setLightWeightPopupEnabled(false);
    }

	public StsComboBoxFieldBean copy(Object beanObject)
	{
		try
		{
			StsComboBoxFieldBean beanCopy = new StsComboBoxFieldBean();
            String labelStg = "";
            if(label != null)
                labelStg = label.getText();
			if(listFieldName != null)
				beanCopy.initialize(beanObject, fieldName, labelStg, listFieldName);
			else
				beanCopy.initialize(beanObject, fieldName, labelStg, items);
			return beanCopy;
		}
		catch(Exception e)
		{
			return null;
		}
	}

    protected void initializeGetItems(Class c, String listFieldName)
    {
		try
		{
			getItems = StsToolkit.getAccessor(c, listFieldName, "get", null);
            if(getItems == null)
            {
                StsException.systemError(this, "initializeGetItems", "Failed to find getter. Need method with name: get" + listFieldName);
                return;
            }
			this.listFieldName = listFieldName;
            comboBox.addPopupMenuListener(this);
        }
		catch (Exception e)
		{
			e.printStackTrace();
		}
    }

    public void setPopupMenuListener(PopupMenuListener listener)
	{
		comboBox.addPopupMenuListener(listener);
	}

    /**  When the view dropdown is selected, first update the list. */
     public void popupMenuWillBecomeVisible(PopupMenuEvent e)
     {
         if(beanObject != null) setItemsFromBeanObject();
//         StsException.systemDebug(this, "popupMenuWillBecomeVisible", "called");
     }

     protected Object[] getItemsFromComboBox()
     {
         ComboBoxModel model = comboBox.getModel();
         int nItems = model.getSize();
         Object[] items = new Object[nItems];
         for(int n = 0; n < nItems; n++)
         {
             items[n] = model.getElementAt(n);
         }
         return items;
     }

     public void setItemsFromBeanObject()
     {
         if(beanObject == null) return;

         try
         {
            if(getItems != null)
            {
                Object[] oldItems = getItemsFromComboBox();
                Object selectedItem = getValueObject();
                Object[] items = (Object[])getItems.invoke(beanObject, new Object[0]);
                if(listsCompare(items, oldItems)) return;
                setListItems(items, selectedItem);
    //          StsException.systemDebug(this, "setItemsFromBeanObject", "number of items: " + items.length);
             }
         }
         catch(Exception ex)
         {
             StsException.outputWarningException(this, "setItemsFromBeanObject()", ex);
         }
     }

     private boolean listsCompare(Object[] newList, Object[] oldList)
     {
         if(newList == null || oldList == null) return false;
         if(newList.length != oldList.length) return false;
         for(int n = 0; n < newList.length; n++)
            if(newList[n] != oldList[n]) return false;
         return true;
     }

    /**  When the view dropdown is selected, first update the list. */
     public void popupMenuWillBecomeInvisible(PopupMenuEvent e)
     {
     }

     /**  When the view dropdown is selected, first update the list. */
     public void popupMenuCanceled(PopupMenuEvent e)
     {

     }

	public void setLightWeightPopupEnabled(boolean enabled)
    {
        comboBox.setLightWeightPopupEnabled(enabled);
    }

    protected void layoutBean()
    {
        try
        {
            setLayout(new BorderLayout());
            add(comboBox, BorderLayout.CENTER);
            if (label != null) add(label, BorderLayout.WEST);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    public void addActionListener()
    {
        comboBox.addActionListener(this);
		comboBox.addPopupMenuListener(this);
//        comboBox.addPopupMenuListener(this);
    }

    public void removeActionListener()
    {
		System.out.println("Remove ");
        comboBox.removeActionListener(this);
        comboBox.removePopupMenuListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
		//System.out.println("action"+ e);
        if (comboBox.getModel().getSize() == 0)return;
        setValueInPanelObjects();
    }

	public void itemStateChanged(ItemEvent e)
	   {
		   System.out.println("item"+ e);
		   //if (comboBox.getModel().getSize() == 0)return;
		   //setValueInPanelObjects();
    }
    public Component[] getBeanComponents()
    {
        return new Component[]
            {label, comboBox};
    }

    public Object getValueObject()
    {
        return comboBox.getSelectedItem();
    }

	public String toString()
	{
        Object selectedItem = comboBox.getSelectedItem();
        if(selectedItem != null)
            return selectedItem.toString();
        else
            return StsParameters.NONE_STRING;
    }

    public Object fromString(String string)
	{
        Object[] items = getItemsFromComboBox();
        if(items == null) return null;
		for(int n = 0; n < items.length; n++)
			if(items[n].toString().equalsIgnoreCase(string)) return items[n];
		return null;
	}

	public String getLabelString()
	{
        Object selectedItem = comboBox.getSelectedItem();
        String selectedItemString = "selected: none";
        String labelString = "no label";
        if(selectedItem != null)
            selectedItemString = selectedItem.toString();
        if(label != null)
            labelString = label.getText();
        return  labelString + " " + selectedItemString;
    }

/*
	public Object getDefaultValueObject()
	{
		return comboBox.getItemAt(0);
    }
*/
//    public Object getValueObject() { return getString(); }
//    public String getString() { return (String) comboBox.getSelectedItem(); }


    public void setListItems(StsObjectRefList objectRefList)
    {
        setListItems(objectRefList.getArrayList(), null);
    }

    public void setListItems(Object[] items)
    {
        setListItems(items, null);
    }

    public void setListItems(Object[] items_, Object selectedItem_)
    {
        // if (items_ == null || items_.length == 0) return;
		setEditable(true);
        final Object[] finalItems = items_;
		this.items = items_;
        final Object selectedItem = selectedItem_;
        final boolean listHasItem = StsMath.arrayContains(finalItems, selectedItem);
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run()
        {
            DefaultComboBoxModel comboModel;

            if(finalItems == null || finalItems.length == 0)
            {
                comboModel = new DefaultComboBoxModel(  );
                comboBox.setModel(comboModel);
                return;
            }

            comboModel = new DefaultComboBoxModel(finalItems);
            comboBox.setModel(comboModel);
            if(selectedItem != null && listHasItem)
                comboBox.setSelectedItem(selectedItem);
            else
                comboBox.setSelectedIndex(0);
            repaint();
        } } );
    }

    public void setSelectedItem(final Object item)
    {
        Runnable runnable = new Runnable() { public void run() { comboBox.setSelectedItem(item); }};
        StsToolkit.runWaitOnEventThread(runnable);
    }

    // override in subclass if object is inside an item
    public Object getItem(Object object)
    {
        //int nItems = getNItems();
        //if (nItems == 0)return null;
        //else return object;
		return object;
    }

    public int getNItems()
    {
        return comboBox.getModel().getSize();
    }

    public void setToolTipText(String tip)
    {
        comboBox.setToolTipText(tip);
    }

    public int getSelectedIndex()
    {
        return comboBox.getSelectedIndex();
    }

    public void setSelectedIndex(int index_)
    {
        final int index = index_;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { comboBox.setSelectedIndex(index); } } );
    }

    public void setToLastItem()
    {
        int nItems = comboBox.getModel().getSize();
        setSelectedIndex(nItems - 1);
    }

    public JComboBox getComboBox()
    {
        return comboBox;
    }

    public void setPreferredSize(int width, int height)
    {
        final int w = width;
        final int h = height;
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.setPreferredSize(new Dimension(w, h)); } } );
    }

    public void addItem(Object object)
    {
        final Object obj = object;
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.setEnabled(editable); comboBox.addItem(obj);} } );
    }

    protected void setEditable()
    {
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.setEnabled(editable); } } );
    }

    public void removeAll()
    {
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.removeAllItems(); } } );
    }

    public void repaint()
    {
        if(comboBox == null) return;
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.repaint(); } } );
    }

    public Component getMainComponent() { return comboBox; }

    public void doSetValueObject(Object object)
    {
        final Object item = getItem(object);
        if (item == null)return;
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.getModel().setSelectedItem(item); } } );
    }

    public Object getValue()
    {
        return comboBox.getSelectedItem();
    }
	public void setValueFromPanelObject(Object panelObject)
	{
        if(panelObject == null) return;
        super.setValueFromPanelObject(panelObject);
        setItemsFromBeanObject();
    }

    public Object getCellEditorValue()
    {
        return comboBox.getSelectedItem();
    }

    public boolean shouldSelectCell(EventObject anEvent)
    {
        if (anEvent instanceof MouseEvent)
        {
            MouseEvent e = (MouseEvent)anEvent;
            return e.getID() != MouseEvent.MOUSE_DRAGGED;
        }
        return editable;
    }

    public boolean stopCellEditing()
    {
		if (comboBox.isEditable())
        {
		    comboBox.actionPerformed(new ActionEvent(this, 0, ""));
		}
		return super.stopCellEditing();
	}
    public void addKeySelector()
    {
       comboBox.setEditable(false);
        // Install the custom key selection manager
       comboBox.setKeySelectionManager(new KeySelectionManager());
    }

    static public void main(String[] args)
    {
        ComboBoxKeySelectionTest test = new ComboBoxKeySelectionTest();
        StsComboBoxFieldBean comboBean = new StsComboBoxFieldBean(test, "selected", "Selection", ComboBoxKeySelectionTest.items);
        comboBean.addKeySelector();
        StsJPanel panel = new StsJPanel();
        panel.add(comboBean);
        StsToolkit.createDialog(panel);
    }
}
/*
    static public void main(String[] args)
    {
        ComboBoxTestObject testObject1 = new ComboBoxTestObject("testObject1");
        ComboBoxTestObject testObject2 = new ComboBoxTestObject("testObject2");
        ComboBoxTestObject[] testObjects = new ComboBoxTestObject[] { testObject1, testObject2 };
        ComboBoxTestPanelObject panelObject = new ComboBoxTestPanelObject();
        StsComboBoxFieldBean comboBoxFieldBean = new StsComboBoxFieldBean(panelObject, "testObject1", "Test Object", testObjects);
		comboBoxFieldBean.setSelectedItem(testObject2);
		StsJPanel panel = new StsJPanel();
        panel.add(comboBoxFieldBean);
        StsToolkit.createDialog(panel);
    }
*/
// This key selection manager will handle selections based on multiple keys.
class KeySelectionManager implements JComboBox.KeySelectionManager
{
    long lastKeyTime = 0;
    String pattern = "";

    public int selectionForKey(char aKey, ComboBoxModel model)
    {
        // Find index of selected item
        int selIx = 01;
        Object sel = model.getSelectedItem();
        if(sel != null)
        {
            for(int i = 0; i < model.getSize(); i++)
            {
                if(sel.equals(model.getElementAt(i)))
                {
                    selIx = i;
                    break;
                }
            }
        }

        // Get the current time
        long curTime = System.currentTimeMillis();

        // If last key was typed less than 300 ms ago, append to current pattern
        if(curTime - lastKeyTime < 300)
        {
            pattern += ("" + aKey).toLowerCase();
        }
        else
        {
            pattern = ("" + aKey).toLowerCase();
        }

        // Save current time
        lastKeyTime = curTime;

        // Search forward from current selection
        for(int i = selIx + 1; i < model.getSize(); i++)
        {
            String s = model.getElementAt(i).toString().toLowerCase();
            if(s.startsWith(pattern))
            {
                return i;
            }
        }

        // Search from top to current selection
        for(int i = 0; i < selIx; i++)
        {
            if(model.getElementAt(i) != null)
            {
                String s = model.getElementAt(i).toString().toLowerCase();
                if(s.startsWith(pattern))
                {
                    return i;
                }
            }
        }
        return -1;
    }
}

class ComboBoxKeySelectionTest
{
    String selected = items[0];

    public static String[] items = {"Ant", "Ape", "Bat", "Boa", "Cat", "Cow"};

    ComboBoxKeySelectionTest()
    {
    }

    public void setSelected(String item)
    {
        selected = item;
        System.out.println("Selected: " + item);
    }

    public String getSelected() { return selected; }
}

class ComboBoxTestPanelObject
{
    ComboBoxTestObject testObject1 = null;
	ComboBoxTestObject testObject2 = null;

    ComboBoxTestPanelObject()
    {
    }

    public void setTestObject1(Object testObject)
    {
        this.testObject1 = (ComboBoxTestObject)testObject;
        System.out.println("setTestObject1 " + this.testObject1.toString());
    }
    public Object getTestObject1() { return testObject1; }

	public void setTestObject2(Object testObject)
	{
		this.testObject2 = (ComboBoxTestObject)testObject;
		System.out.println("setTestObject2 " + this.testObject2.toString());
	}
    public Object getTestObject2() { return testObject2; }
}

class ComboBoxTestObject extends StsMainObject
{
    ComboBoxTestObject(String name)
    {
        super(false);
        setName(name);
    }

//    public String toString() { return name; }
}
