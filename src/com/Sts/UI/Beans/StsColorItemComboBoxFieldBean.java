//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Beans;

import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/** This class is intended for a list of instances of StsColorListItem.
 *  Each instance has an associated object, name, and color.
 */
public class StsColorItemComboBoxFieldBean extends StsComboBoxFieldBean
{
//	private StsColorListItem[] items = null;

	public StsColorItemComboBoxFieldBean()
    {
    }

    public StsColorItemComboBoxFieldBean(String fieldLabel)
    {
        classInitialize(null, null, fieldLabel);
        comboBox.setRenderer(new ColorItemListRenderer());
//        comboBox.setLightWeightPopupEnabled(false);
    }

    public StsColorItemComboBoxFieldBean(String fieldName, String fieldLabel)
    {
        classInitialize(null, fieldName, fieldLabel);
        comboBox.setRenderer(new ColorItemListRenderer());
//        comboBox.setLightWeightPopupEnabled(false);
    }

	public StsColorItemComboBoxFieldBean(Class c, String fieldName, String fieldLabel, StsColorListItem[] items)
    {
//        this.items = items;
        initialize(c, fieldName, fieldLabel, items);
        comboBox.setRenderer(new ColorItemListRenderer());
//        comboBox.setLightWeightPopupEnabled(false);
   }

	public StsColorItemComboBoxFieldBean(Class c, String fieldName, String fieldLabel)
    {
        super.classInitialize(c, fieldName, fieldLabel, null);
        comboBox.setRenderer(new ColorItemListRenderer());
//        comboBox.setLightWeightPopupEnabled(false);
   }

	public StsColorItemComboBoxFieldBean(Object beanObject, String fieldName, String fieldLabel, StsColorListItem[] items)
    {
        initialize(beanObject, fieldName, fieldLabel, items);
        comboBox.setRenderer(new ColorItemListRenderer());
//        comboBox.setLightWeightPopupEnabled(false);
   }

    public void initialize(Class c, String fieldName, String fieldLabel, StsColorListItem[] items)
    {
        super.classInitialize(c, fieldName, fieldLabel, items);
        comboBox.setRenderer(new ColorItemListRenderer());
    }

    public void initialize(Object beanObject, String fieldName, String fieldLabel, StsColorListItem[] items)
    {
        if(beanObject == null)
        {
            StsException.systemError("StsFieldBean.classInitialize() failed. beanObject cannot be null.");
            return;
        }
        this.beanObject = beanObject;
        Class c = beanObject.getClass();
        initialize(c, fieldName, fieldLabel, items);
        setValueFromPanelObject(beanObject);
    }

    public void initialize(Object instance, String fieldName, String fieldLabel)
    {
        super.initialize(instance, fieldName, fieldLabel);
        comboBox.setRenderer(new ColorItemListRenderer());
    }

	public StsColorItemComboBoxFieldBean copy(Object beanObject)
	{
		try
		{
			StsColorItemComboBoxFieldBean beanCopyItem = new StsColorItemComboBoxFieldBean();
			if(listFieldName != null)
				beanCopyItem.initialize(beanObject, fieldName, getLabelString(), listFieldName);
			else
				beanCopyItem.initialize(beanObject, fieldName, getLabelString(), items);
			return beanCopyItem;
		}
		catch(Exception e)
		{
			return null;
		}
	}
    
    public void initializeColors(String fieldLabel, StsColorListItem[] items)
    {
        super.classInitialize(null, null, fieldLabel, items);
        comboBox.setRenderer(new ColorItemListRenderer());
    }

    public void setListItems(StsColorListItem[] items)
    {
    	super.setListItems(items);
    }

    public void addItem(StsColorListItem item)
    {
        super.addItem(item);
    }

    public Object getValueObject() { return comboBox.getSelectedItem(); }
//	public Object getDefaultValueObject() { return comboBox.getItemAt(0); }

    public void setValueObject(StsColorListItem item)
    {
        super.setSelectedItem(item);
    }

    public void setSelectedIndex(int index)
    {
        super.setSelectedIndex(index);
    }

	public String toString()
	{
		return comboBox.getSelectedItem().toString();
	}

	public Object fromString(String string)
	{
        if(getNItems() == 0) return null;
		return Color.decode(string);
	}

    public Object getSelectedItem() { return comboBox.getSelectedItem(); }
    public void setSelectedItem(Object object)
    {
        final Object obj = object;
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.setSelectedItem(obj); } } );
    }

    public void removeItem(Object object)
    {
        final Object obj = object;
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.removeItem(obj); } } );
    }

    public void removeItemAtIndex(int index_)
    {
       final int index = index_;
       StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.removeItemAt(index); } } );
    }

    public void removeAllItems()
    {
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { comboBox.removeAllItems(); } } );
    }

    public int getListSize()
    {
        return comboBox.getModel().getSize();
    }

    public StsColorListItem getItemAt(int index)
    {
        return (StsColorListItem)comboBox.getModel().getElementAt(index);
    }

    public void deleteObject(Object object)
    {
        int nItems = getListSize();
        for(int n = 0; n < nItems; n++)
        {
            StsColorListItem item = getItemAt(n);
            if(item.getObject() == object)
            {
                removeItem(item);
                return;
            }
        }
    }

    static public StsColorListItem getNullListItem()
    {
        StsColor color = StsColor.GRAY;
        String name = "";
        return new StsColorListItem(null, color, name, 16, 16, 1);
    }

    static public void main(String[] args)
    {
        StsColorItemListTest test = new StsColorItemListTest();
        StsColorItemComboBoxFieldBean colorItemListBean = new StsColorItemComboBoxFieldBean();
        colorItemListBean.initialize(test, "colorListItem", "Colors", new Object[] { StsColorListItem.nullColorListItem() } );
        test.initialize(colorItemListBean);
		StsButton addColorButton = new StsButton("addColor", "Add a color.", test, "addColor", colorItemListBean);
        StsButton deleteAllButton = new StsButton("deleteAll", "delete all colors.", test, "deleteAllColors", colorItemListBean);
		StsJPanel panel = StsJPanel.addInsets();
		panel.add(colorItemListBean);
 		panel.add(addColorButton);
        panel.add(deleteAllButton);
        StsToolkit.createDialog(panel);
    }
}

class ColorItemListRenderer extends JLabel implements ListCellRenderer
{
    public Component getListCellRendererComponent(JList list, Object value, int index,
                        boolean isSelected, boolean hasFocus)
    {
        StsColorListItem item;

        if(value == null)
            item = StsColorItemComboBoxFieldBean.getNullListItem();
        else
            item = (StsColorListItem)value;

        setOpaque(true);
        setIcon(item.getIcon());
        setText(item.getName());
        if (isSelected)
        {
            setBackground(Color.black);
            setForeground(Color.white);
        }
        else // not selected
        {
            setBackground(Color.white);
            setForeground(Color.black);
        }
        return this;
    }
}

class StsColorItemListTest
{
    public StsColorListItem colorListItem = StsColorListItem.nullColorListItem();
    public int index = 0;

    public static StsColor[] colors = StsColor.colors8;

    StsColorItemListTest()
    {
    }

//    public Color[] getColors() { return colors; }
//    public Color getColor() { return color; }
//    public void setBeachballColors(Color c) { color = c; }

    public StsColorListItem getColorListItem() { return colorListItem; }
    public void setColorListItem(StsColorListItem colorListItem) { this.colorListItem = colorListItem; }

    public void initialize(StsColorItemComboBoxFieldBean colorItemListBean)
    {
        StsColorListItem listItem = StsColorItemComboBoxFieldBean.getNullListItem();
        colorItemListBean.addItem(listItem);
        colorItemListBean.comboBox.setSelectedItem(listItem);
    }

    public void addColor(StsColorItemComboBoxFieldBean colorItemListBean)
    {

        StsColor color = colors[(index++)%10];
        String name = Integer.toString(index-1);

        // first item is a dummy and has a null object; if this is only item,
        // change it to a legitimate item
        StsColorListItem listItem = (StsColorListItem) colorItemListBean.getSelectedItem();
        if(listItem.getObject() == null)
        {
            listItem.setStsColor(color);
            listItem.setName(name);
            listItem.setObject(new String(""));
            colorItemListBean.repaint();
        }
        else
        {
            listItem = new StsColorListItem(new String(""), color, name, 16, 16, 1);
            colorItemListBean.addItem(listItem);
        }
        colorItemListBean.comboBox.setSelectedItem(listItem);
    }

    public void deleteAllColors(StsColorItemComboBoxFieldBean colorItemListBean)
    {
        colorItemListBean.removeAllItems();
        initialize(colorItemListBean);
        colorItemListBean.repaint();
    }
}

