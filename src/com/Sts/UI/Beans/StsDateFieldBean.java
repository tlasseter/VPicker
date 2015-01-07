package com.Sts.UI.Beans;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.util.*;

public class StsDateFieldBean extends StsTextFieldBean implements MouseListener, ChangeListener
{
    String time = "";
    Date date = new Date();
    static SimpleDateFormat format = new SimpleDateFormat("dd-MM-yy HH:mm:ss.S");
    JPopupMenu popup = null;

    private static final int NONE = -1;
    private static final int DAY = 2;
    private static final int MONTH = 5;
    private static final int YEAR = 8;
    private static final int HOUR = 11;
    private static final int MINUTE = 14;
    private static final int SECOND = 17;
    private int currentField = NONE;

    public StsDateFieldBean()
    {
    }

    public StsDateFieldBean(String label)
    {
        this(true, label);
    }

    public StsDateFieldBean(boolean editable)
    {
        this(editable, null);
    }

    public StsDateFieldBean(boolean editable, String label)
    {
        this.classInitialize((Class)null, null, editable, label);
    }

    public StsDateFieldBean(String fieldName, String label)
    {
        this.initialize(null, fieldName, time, true, label);
    }

    public StsDateFieldBean(Class c, String fieldName)
    {
        this.initialize(c, fieldName, time, true, null);
    }

    public StsDateFieldBean(Class c, String fieldName, String label)
    {
        this.initialize(c, fieldName, time, true, label);
    }

    public StsDateFieldBean(Class c, String fieldName, boolean editable, String label)
    {
        this.initialize(c, fieldName, time, editable, label);
    }

    public StsDateFieldBean(Class c, String fieldName, boolean editable)
    {
        this.initialize(c, fieldName, time, editable, "");
    }

    public StsDateFieldBean(Class c, String fieldName, String time, String label)
    {
        this.initialize(c, fieldName, time, true, label);
    }

	public StsDateFieldBean(Class c, String fieldName, String time, boolean editable, String label)
	{
		this.initialize(c, fieldName, time, editable, label);
	}

	public StsDateFieldBean(Object beanObject, String fieldName, String label)
	{
		this(beanObject, fieldName, true, label);
		layoutBean();
    }

    public StsDateFieldBean(Object beanObject, String fieldName, boolean editable, String label)
    {
        initialize(beanObject, fieldName, editable, label);
		layoutBean();
    }

    public StsDateFieldBean(Object beanObject, String fieldName, String value, boolean editable, String label)
    {
        this.initialize(beanObject, fieldName, value, editable, label);
    }

    public void initialize(Class c, String fieldName, String value, String label)
    {
        super.classInitialize(c, fieldName, true, label);
		initializeValue(value);
		layoutBean();
   }

   private void initializeValue(String value)
   {
	   this.time = value;
	   super.setText(value);
   }

    public void initialize(Class c, String fieldName, String value, boolean editable, String label)
    {
        super.classInitialize(c, fieldName, editable, label);
        initializeValue(value);
		layoutBean();
    }

	public void classInitialize(Class c, String fieldName, boolean editable, String label)
	{
		super.classInitialize(c, fieldName, editable, label);
		layoutBean();
	}

    public void initialize(Object beanObject, String fieldName, String value, boolean editable, String label)
    {
        super.initialize(beanObject, fieldName, editable, label);
        initializeValue(value);
		layoutBean();
    }

	public void initialize(Object beanObject, String fieldName, boolean editable, String label)
	{
		super.initialize(beanObject, fieldName, editable, label);
		layoutBean();
    }

	public StsDateFieldBean copy(Object beanObject)
	{
		StsDateFieldBean beanCopy = new StsDateFieldBean();
		beanCopy.initialize(beanObject, fieldName, editable, getLabelText());
		return beanCopy;
	}

    public void setFormat(SimpleDateFormat newFormat)
    {
        format = newFormat;
    }

	protected void layoutBean()
	{
		textField.setColumns(10);
        textField.setBorder(BorderFactory.createLineBorder(Color.GRAY));
        textField.setForeground(Color.GREEN);
        textField.setBackground(Color.BLACK);
        textField.setFont(new Font("Dialog", Font.BOLD, 11));
        textField.setToolTipText("Enter model date and time.");
        if(editable)
            textField.addMouseListener(this);
		labelAndComponentLayout(label, textField);
    }

    // return true if value has changed
    public boolean setValue(String newValue)
    {
        try
        {
            String oldValue = time;
            time = newValue;
            setText(newValue);
            return !newValue.equals(oldValue);
        }
        catch(Exception e)
        {
            StsToolkit.beep();
            textField.requestFocus();
            return false;
        }
    }

	public String toString() { return time; }
	public Object fromString(String string) { return string; }

    public boolean setValueObject(Object object)
    {
        if(!(object instanceof String)) return false;
        return setValue((String)object);
    }

    public Object getValueObject()
    {
        return time;
    }

    public String getValue()
    {
        return textField.getText();
    }

    protected boolean isStringBean() { return true; }

    public void mouseClicked(MouseEvent e)
    {
        int clicks = e.getClickCount();
        if(clicks == 2)
        {
            int selectedEnd = textField.getSelectionEnd();
            String selectedValue = textField.getSelectedText();
            if(selectedValue.equals(":") || selectedValue.equals("-"))
                return;
            if(selectedValue.startsWith("0"))
                selectedValue = selectedValue.substring(1,2);
            if(selectedValue.length() > 2)
                selectedValue = selectedValue.substring(0,2);

            int value = -1;
            try { value = Integer.parseInt(selectedValue); }
            catch(Exception ex) { StsMessageFiles.infoMessage("Invalid Date Input."); return; }
            if(value == -1)
               return;

            String subString = textField.getText().substring(0,selectedEnd);
            switch(subString.length())
            {
                case DAY:
                    currentField = DAY;
                    showSlider("Day", e.getX(), e.getY(), 1, 31, 5, 0, value);
                    break;
                case MONTH:
                    currentField = MONTH;
                    showSlider("Month", e.getX(), e.getY(), 1, 12, 1, 0, value);
                    break;
                case YEAR:
                    currentField = YEAR;
                    Calendar c = Calendar.getInstance();
                    int year = c.get(c.YEAR);
                    if(value > 25)
                        value = value + 1900;
                    else
                        value = value + 2000;
                    showSlider("Year", e.getX(), e.getY(), 1969, year, 10, 5, value);
                    break;
                case HOUR:
                    currentField = HOUR;
                    showSlider("Hour", e.getX(), e.getY(), 0, 23, 4, 2, value);
                    break;
                case MINUTE:
                    currentField = MINUTE;
                    showSlider("Minute", e.getX(), e.getY(), 0, 60, 10, 5, value);
                    break;
                default:
                    currentField = SECOND;
                    showSlider("Second", e.getX(), e.getY(), 0, 60, 10, 5, value);
                    break;
            }
        }
    }
    public void mousePressed(MouseEvent e) { }
    public void mouseReleased(MouseEvent e)
    {
        if(e.getSource() instanceof JSlider)
        {
            popup.setVisible(false);
        }
    }
    public void mouseEntered(MouseEvent e) { }
    public void mouseExited(MouseEvent e) { }

    private void showSlider(String label, int xpos, int ypos, int minValue, int maxValue, int major, int minor, int value)
    {
        popup = new JPopupMenu();
        popup.setLayout(new GridBagLayout());
        JLabel jlabel = new JLabel(label);
        jlabel.setFont(new Font("Dialog", Font.BOLD, 10));
        popup.add(jlabel);
        popup.add(jlabel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0 ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

        JSlider slider = new JSlider(JSlider.VERTICAL, minValue, maxValue, value);
        slider.setInverted(true);
        slider.setPaintTrack(false);
        slider.setPaintLabels(true);
        slider.setPaintTicks(true);
        slider.setMajorTickSpacing(major);
        slider.setMinorTickSpacing(minor);
        slider.setFont(new Font("Dialog", Font.PLAIN, 10));
        slider.setBackground(Color.LIGHT_GRAY);
        slider.setBorder(BorderFactory.createEtchedBorder());
        slider.addChangeListener(this);
        slider.addMouseListener(this);
        popup.add(slider,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0 ,GridBagConstraints.WEST, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));

        textField.add(popup);
        popup.show(textField, xpos + 12, ypos);
    }

    public void stateChanged(ChangeEvent e)
    {
        if(e.getSource() instanceof JSlider)
        {
            String sValue = null;
            int value = ((JSlider)e.getSource()).getValue();
            String cTime = textField.getText();
            sValue = Integer.toString(value);
            if(sValue.length() == 1)
                sValue = "0" + sValue;
            if(currentField == YEAR)
                sValue = sValue.substring(2,4);
            cTime = cTime.substring(0,currentField-2) + sValue + cTime.substring(currentField,cTime.length());

            try
            {
                Date date = format.parse(cTime);
                if (date.getTime() >= System.currentTimeMillis())
                {
                    date.setTime(System.currentTimeMillis());
                    time = format.format(date);
                }
                else
                {
                    textField.setText(cTime);
                    time = cTime;
                }
                setValueInPanelObject();
            }
            catch(Exception ex)
            {
                StsMessageFiles.infoMessage("Unable to parse date: " + time);
            }
        }
    }

    static public String convertToString(long time)
    {
        return format.format(new Date(time));
    }

    static public long convertToLong(String sDate)
    {
        try
        {
            Date date = format.parse(sDate);
            return date.getTime();
        }
        catch(Exception e)
        {
            StsMessageFiles.infoMessage("Unable to parse date: " + sDate);
            return -1l;
        }
    }
    static public boolean validateDateInput(String sDate)
    {
        try
        {
            Date date = format.parse(sDate);
            return true;
        }
        catch(Exception e)
        {
            return false;
        }
    }
}
