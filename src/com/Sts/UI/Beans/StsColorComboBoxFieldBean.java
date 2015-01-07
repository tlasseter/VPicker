package com.Sts.UI.Beans;

import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsColorComboBoxFieldBean extends StsComboBoxFieldBean
{
	private StsColorListItem[] items = null;

	public StsColorComboBoxFieldBean()
    {
    }

    public StsColorComboBoxFieldBean(String fieldLabel)
    {
        super.initialize(null, null, true, fieldLabel);
    }

    public StsColorComboBoxFieldBean(String fieldName, String fieldLabel)
    {
        super.initialize(null, fieldName, true, fieldLabel);
    }

	public StsColorComboBoxFieldBean(Class c, String fieldName, String fieldLabel, StsColor[] stsColors)
    {
        initializeColors(c, fieldName, fieldLabel, stsColors);
    }

	public StsColorComboBoxFieldBean(Class c, String fieldName, String fieldLabel, StsSpectrum spectrum)
    {
        initializeColors(c, fieldName, fieldLabel, spectrum.getStsColors());
    }

    public StsColorComboBoxFieldBean(Class c, String fieldName, String fieldLabel)
    {
        super.classInitialize(c, fieldName, fieldLabel, null);
    }

	public StsColorComboBoxFieldBean(Object beanObject, String fieldName, String fieldLabel, StsColor[] colors)
    {
        initializeColors(beanObject, fieldName, fieldLabel, colors);
    }

	public StsColorComboBoxFieldBean copy(Object beanObject)
	{
		try
		{
			StsColor[] colors = StsColorListItem.getColorsFromItems(items);
			return new StsColorComboBoxFieldBean(beanObject, fieldName, getLabelText(), colors);
		}
		catch(Exception e)
		{
			return null;
		}
	}

    public void initializeColors(Class c, String fieldName, String fieldLabel, StsSpectrum spectrum)
    {
        StsColor[] colors = spectrum.getStsColors();
        super.classInitialize(c, fieldName, fieldLabel, createItems(colors));
        comboBox.setRenderer(new StsColorListRenderer());
    }

    public void initializeColors(Object beanObject, String fieldName, String fieldLabel, StsSpectrum spectrum)
    {
        StsColor[] colors = spectrum.getStsColors();
        if(beanObject == null)
        {
            StsException.systemError("StsFieldBean.classInitialize() failed. beanObject cannot be null.");
            return;
        }
        this.beanObject = beanObject;
        Class c = beanObject.getClass();
        initializeColors(c, fieldName, fieldLabel, colors);
        setValueFromPanelObject(beanObject);
    }

    public void initializeColors(Class c, String fieldName, String fieldLabel, StsColor[] colors)
    {
        super.classInitialize(c, fieldName, fieldLabel, createItems(colors));
        comboBox.setRenderer(new StsColorListRenderer());
    }

	public void initialize(Object instance, String fieldName, String fieldLabel, StsColorListItem[] colorListItems)
    {
        super.initialize(instance, fieldName, fieldLabel, colorListItems);
        comboBox.setRenderer(new StsColorListRenderer());
    }


    public void initializeColors(Object beanObject, String fieldName, String fieldLabel, StsColor[] colors)
    {
        if(beanObject == null)
        {
            StsException.systemError("StsFieldBean.classInitialize() failed. beanObject cannot be null.");
            return;
        }
        this.beanObject = beanObject;
        Class c = beanObject.getClass();
        initializeColors(c, fieldName, fieldLabel, colors);
        setValueFromPanelObject(beanObject);
    }

    public void initializeColors(String fieldLabel, StsColor[] colors)
    {
        super.classInitialize(null, null, fieldLabel, createItems(colors));
        comboBox.setRenderer(new StsColorListRenderer());
    }

    private StsColorListItem[] createItems(StsColor[] colors)
    {
    	if( colors == null ) return null;
        items = new StsColorListItem[colors.length];
        for( int i=0; i<items.length; i++ )
        {
			items[i] = new StsColorListItem(colors[i], null, 48, 16, 1);
        }
        return items;
    }

    public void setListItems(StsColor[] colors)
    {
    	super.setListItems(createItems(colors));
    }

    public void setListItems(StsSpectrum spectrum)
    {
        setListItems(spectrum.getStsColors());
    }

    public void add(StsColor color)
    {
        StsColorListItem colorListItem = new StsColorListItem(color, null, 48, 16, 1);
        comboBox.addItem(colorListItem);
    }

    public Object getValueObject() { return getStsColor(); }
/*
	public Object getDefaultValueObject()
	{
		Object listItem = comboBox.getItemAt(0);
		if(listItem == null) return null;
		if(listItem instanceof Color) return (Color)listItem;
		StsColorListItem item = (StsColorListItem)listItem;
		return item == null ? null : item.getColor();
	}
*/
    public StsColor getStsColor()
    {
        Object listItem = comboBox.getSelectedItem();
        if(listItem instanceof StsColor) return (StsColor)listItem;
        StsColorListItem item = (StsColorListItem)listItem;
    	return item == null ? null : item.getStsColor();
	}

    public void setStsColor(StsColor color)
    {
        this.comboBox.setSelectedItem(color);
    }

    public void setValueObject(StsColor color)
    {
        removeActionListener();
    	if( items == null )
        {
            setListItems(new StsColor[] { color });
        }

        for( int i=0; i<items.length; i++ )
        {
			if( items[i].getStsColor().equals(color) )
		    	comboBox.setSelectedIndex(i);
        }
        addActionListener();
        comboBox.repaint();
    }

    public Object getItem(Object object)
    {
        if(items == null) return null;
		if(!(object instanceof StsColor)) return null;
        for(int n = 0; n < items.length; n++)
        {
            StsColorListItem item = items[n];
            StsColor itemColor = item.getStsColor();
			if(itemColor.equals((StsColor)object)) return item;
//            if(colorsEqual(itemColor, objectColor)) return item; // double-check, color-by-color
//            if(colorItem.getColor().equals((Color)object)) return colorItem;
        }
        return null;
    }

    private boolean colorsEqual(Color color0, Color color1)
    {
        float[] rgb0 = new float[4];
        color0.getComponents(rgb0);
        float[] rgb1 = new float[4];
        color1.getComponents(rgb1);
        for(int n = 0; n < 4; n++)
            if(rgb0[n] != rgb1[n]) return false;
        return true;
    }

/*
    public void setValue(Object obj)
    {
    	try
        {
	        Color color = (Color) get.invoke(obj, new Object[0]);
    		setValue(color);
        }
        catch(Exception e) { e.printStackTrace(); }
    }

    public int getInt()
    {
    	return list.getSelectedIndex();
    }
*/


    public String toString()
	{
        Object item = comboBox.getSelectedItem();
        if(item == null)
            return StsParameters.NONE_STRING;
        else
            return item.toString();
	}

	public Object fromString(String string)
	{
        StsColor color = StsColor.colorFromString(string);
		for (int n = 0; n < items.length; n++)
		{
			StsColor itemColor = items[n].getStsColor();
			if(color.equals(itemColor))
				return itemColor;
		}
		return StsColor.BLACK;
	}

    static public void main(String[] args)
    {
        StsColorListTestNew testObject = new StsColorListTestNew();

        StsColor[] colors = testObject.getColors();

		StsColorComboBoxFieldBean colorComboBoxBean1 = new StsColorComboBoxFieldBean();
		colorComboBoxBean1.initializeColors(testObject, "stsColor", "Colors", colors);

        StsColorComboBoxFieldBean colorComboBoxBean = new StsColorComboBoxFieldBean(testObject, "stsColor", "Colors", colors);
        testObject.addColor(colorComboBoxBean);
        testObject.addColor(colorComboBoxBean);
        testObject.addColor(colorComboBoxBean);

        StsButton addColorButton = new StsButton("addColor", "Add a color.", testObject, "addColor", colorComboBoxBean);

		StsJPanel panel = StsJPanel.addInsets();
		panel.add(colorComboBoxBean1);
		panel.add(colorComboBoxBean);
		panel.add(addColorButton);

        StsToolkit.createDialog(panel);
    }
}

class StsColorListTestNew
{
    public StsColor stsColor;
    public int index = 0;

    public static StsColor[] colors = StsColor.colors16;

    StsColorListTestNew()
    {
    }

    public StsColor[] getColors() { return colors; }
    public StsColor getStsColor() { return stsColor; }
    public void setStsColor(StsColor c) { stsColor = c; }

    public void addColor(StsColorComboBoxFieldBean colorComboBoxBean)
    {
        StsColor color = colors[(index++)%16];
        colorComboBoxBean.add(color);
        colorComboBoxBean.comboBox.setSelectedIndex(index-1);
    }

    void printTestValues()
    {
        System.out.println("color: " + stsColor.toLabelString());
    }
}

