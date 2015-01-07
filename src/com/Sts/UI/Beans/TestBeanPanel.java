package com.Sts.UI.Beans;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class TestBeanPanel extends StsFieldBeanPanel
{
	protected Object panelObject = null; // Object from which this panel get/set(s) values
	protected boolean checkbox = false;
	JTextField textField1 = new JTextField(5);
	JTextField textField2 = new JTextField(5);
	JTextField textField3 = new JTextField(5);
	JTextField textField4 = new JTextField(5);
	JTextField textField5 = new JTextField(5);
	StsIntFieldBean intBean1 = new StsIntFieldBean();
	StsFloatFieldBean floatBean2 = new StsFloatFieldBean();
	StsComboBoxFieldBean comboListBean = new StsComboBoxFieldBean();
    StsListFieldBean listBean = new StsListFieldBean();
    StsColorComboBoxFieldBean colorComboBoxBean = new StsColorComboBoxFieldBean();
	StsBooleanFieldBean booleanBean;
	StsRadioButtonFieldBean radioBean = new StsRadioButtonFieldBean();
    ButtonGroup buttonGroup = new ButtonGroup();
    StsStringFieldBean stringBean;
	StsStringFieldBean noEditStringBean = new StsStringFieldBean();
	StsStringFieldBean thisStringBean = new StsStringFieldBean(this, "thisName", "This name");
	StsFloatFieldBean panelFloatBean = new StsFloatFieldBean(this, "panelFloat", "Panel float");

    ListObject listItem;
    float[] numberList2 = new float[2];
    int selectedIndex = -1;

    float panelFloat = 123456789f;

	String thisName = "thisName";

    public TestBeanPanel(Object panelObject)
    {
        this.panelObject = panelObject;

        try
        {
			gbc.anchor = GridBagConstraints.WEST;
			gbc.fill = GridBagConstraints.HORIZONTAL;
			textField1.setHorizontalAlignment(JTextField.TRAILING);
			textField1.setText("1234567890TRAILING");
			textField2.setHorizontalAlignment(JTextField.LEADING);
			textField2.setText("1234567890LEADING");
			textField3.setHorizontalAlignment(JTextField.LEFT);
			textField3.setText("1234567890LEFT");
			textField4.setHorizontalAlignment(JTextField.RIGHT);
			textField4.setText("1234567890RIGHT");
            intBean1.initialize(panelObject, "number1", 0, 10, "number1");
            floatBean2.initialize(this, "number2", true, "number2");
            ListObject[] listItems = new ListObject[] { new ListObject("item1"), new ListObject("item2") };
            comboListBean.initialize(this, "comboListItem", "Combo Items", listItems);
            listBean.initialize(this, "listItem", "List Items", listItems);
            StsColor[] colors = new StsColor[] { StsColor.RED, StsColor.GREEN, StsColor.BLUE};
            colorComboBoxBean.initializeColors(panelObject, "color3", "Colors", colors);
            booleanBean = new StsBooleanFieldBean(this, "checkBox", "Checkbox");
            radioBean.initialize(panelObject, "radio", "Radio Button", buttonGroup);
            radioBean.setText("pushed");
            radioBean.setUnselectedText("pulled");
			stringBean = new StsStringFieldBean(panelObject, "name", "Name");
			noEditStringBean.initialize(panelObject, "noEdit", false, "No Edit");
//			panelFloatBean.classInitialize(this, "panelFloat", "Panel float");

            jbInit();

			floatBean2.setValue(1.e30);
			Thread.sleep(2000);
			floatBean2.setValue(1000);
			Thread.sleep(2000);
			floatBean2.setValue(0.001);
        }
        catch(Exception ex)
		{
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception
    {
	//		gbc.anchor = GridBagConstraints.EAST;
		add(textField1);
		add(textField2);
		add(textField3);
		add(textField4);
		gbc.gridy++;
		add(intBean1);
		intBean1.setMaximumSize(new Dimension(50, 10));
		add(floatBean2);
		add(comboListBean);
        add(listBean);
        add(colorComboBoxBean);
		add(booleanBean);
		add(radioBean);
		add(stringBean);
		add(noEditStringBean);
		add(thisStringBean);
		add(panelFloatBean);
    }

    public Object getPanelObject() { return panelObject; }

    public Object[] getPanelObjects() { return new Object[] { panelObject }; }
//    public Object setPanelObject() { return panelObject; }

    public void add(StsFieldBean fieldBean, Object constraints)
    {
        fieldBean.setBeanPanel(this);
        super.add(fieldBean, constraints);
    }

	public void setThisName(String name) { thisName = name; }
	public String getThisName() { return thisName; }

	public boolean getCheckBox() { return checkbox; }
	public void setCheckBox(boolean b)
	{
		intBean1.setEditable(b);
		floatBean2.setEditable(b);
		checkbox = b;
    }

    public ListObject getListItem() { return listItem; }
    public void setListItem(ListObject item)
    {
        listItem = item;
        selectedIndex = listBean.getSelectedIndex();
        floatBean2.setValue(numberList2[selectedIndex]);
        System.out.println("setListItem called: index = " + selectedIndex + " value set to " + numberList2[selectedIndex]);
    }

    public ListObject getComboListItem() { return listItem; }
    public void setComboListItem(ListObject item)
    {
        listItem = item;
        int index = comboListBean.getSelectedIndex();
        floatBean2.setValue(numberList2[index]);
        System.out.println("setCombotListItem called: index = " + index + " value set to " + numberList2[index]);
    }

    public float getNumber2()
    {
        return numberList2[selectedIndex];
    }
    public void setNumber2(float f)
    {
        numberList2[selectedIndex] = f;
        System.out.println("setNumber2 called: index = " + selectedIndex + " value set to " + numberList2[selectedIndex]);
    }

    static public void main(String[] args)
    {
		TestObject panelObject = new TestObject();
		TestBeanPanel panel = new TestBeanPanel(panelObject);
		panel.setMaximumSize(new Dimension(100, 100));
		StsToolkit.createDialog(panel, false, 100, 100);
		System.out.println("checkBox: " + panel.checkbox);
		panelObject.printTestValues();
    }

	public void setPanelFloat(float value) { this.panelFloat = value; }
	public float getPanelFloat() { return panelFloat; }
}

class ListObject
{
    String name;

    public ListObject(String name)
    {
        this.name = name;
    }

    public String toString() { return name; }
}

class TestObject
{
    int number1 = 0;
    float number2 = 60f;
    ListObject listItem;
    Color color3 = Color.green;
//    boolean checkBox;
    boolean radio;
	String name = "Tommy";
	String noEdit = "No edit";

    TestObject()
    {
    }

    public int getNumber1() { return number1; }
    public void setNumber1(int n) { number1 = n; }
//    public float getNumber2() { return number2; }
//    public void setNumber2(float f) { number2 = f; }
    public Color getColor3() { return color3; }
    public void setColor3(Color c) { color3 = c; }
    public ListObject getListItem() { return listItem; }
    public void setListItem(ListObject item) { listItem = item; }
    public boolean getRadio() { return radio; }
    public void setRadio(boolean b) { radio = b; }
	public String getName() { return name; }
	public void setName(String string) { this.name = string; }
	public String getNoEdit() { return noEdit; }

    void printTestValues()
    {
        System.out.println("number1: " + getNumber1());
        System.out.println("number2: " + number2);
        System.out.println("listItem: " + listItem);
        System.out.println("color3: " + color3.toString());
        System.out.println("radio: " + radio);
		System.out.println("name: " + name);
		System.out.println("noEdit: " + noEdit);
    }
}


