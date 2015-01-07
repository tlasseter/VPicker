package com.Sts.UI.Beans;

import com.Sts.DBTypes.StsColor;
import com.Sts.Utilities.StsToolkit;

import javax.swing.*;
import java.awt.*;

/**
 * Copyright:  Copyright (c) 2011
 * Author: Tom Lasseter
 * Date: 11/7/11
 */
public class StsNameColorComboListTest
{
	static public void main(String[] args)
	{
		TestNameColorCombo testList = new TestNameColorCombo();
		StsJPanel listPanel = new StsJPanel();
		StsJPanel[] panelItems =  testList.panelItems;
		for(StsJPanel panelItem : panelItems)
			listPanel.add(panelItem);
		JScrollPane itemScrollPane = new JScrollPane();
		itemScrollPane.getViewport().add(listPanel, null);
		StsJPanel itemsBox = new StsGroupBox("Box of items");
		itemsBox.setPreferredSize(200, 100);
		itemsBox.gbc.fill = GridBagConstraints.BOTH;
		itemsBox.addEndRow(itemScrollPane);
		StsJPanel panel = new StsJPanel();
		panel.gbc.fill = GridBagConstraints.BOTH;
		panel.add(itemsBox);
		StsToolkit.createDialog(panel);
		testList.print();
	}
}

class TestNameColorCombo
{
	StsJPanel[] panelItems;
	StsColor[] colors = StsColor.colors32;
	String[] colorNames = StsColor.colorNames32;
	NameColorObject[] colorItemObjects = new NameColorObject[32];

	TestNameColorCombo()
	{
		panelItems = new StsJPanel[32];
		for(int n = 0; n < 32; n++)
		{
			NameColorObject colorItemObject = new NameColorObject(colorNames[n], colors[n], colors);
			colorItemObjects[n] = colorItemObject;
			StsJPanel panelItem = new StsJPanel();
			  StsStringFieldBean nameBean = new StsStringFieldBean(colorItemObject, "name", true);
			panelItem.addToRow(nameBean);
			StsColorComboBoxFieldBean colorBean = new StsColorComboBoxFieldBean(colorItemObject, "color", null, colors);
			colorBean.setSelectedIndex(n);
			panelItem.addToRow(colorBean);
			panelItems[n] = panelItem;
		}
	}

	void print()
	{
		for(int n = 0; n < 32; n++)
			System.out.println("item " + n + " name " + colorNames[n] + " color Index " + colorItemObjects[n].colorIndex);
	}
}

class NameColorObject
{
	private String name;
	private StsColor color;
	private StsColor[] colors;
	public int colorIndex;
	StsStringFieldBean nameBean;
	StsComboBoxFieldBean colorBean;

	NameColorObject(String name, StsColor color, StsColor[] colors)
	{
		this.name = name;
		this.color = color;
		this.colors = colors;
		nameBean = new StsStringFieldBean(this, "name", true);
		colorBean = new StsComboBoxFieldBean(this, "color", colors);
	}

	public String getName()
	{
		return name;
	}

	public void setName(String name)
	{
		this.name = name;
	}

	public StsColor getColor()
	{
		return color;
	}

	public void setColor(StsColor color)
	{
		this.color = color;
		this.colorIndex = colorBean.getSelectedIndex();
		System.out.println("colorIndex " + colorIndex + " color " + color.toString());
	}
}
