package com.Sts.UI.Beans;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsLongFieldBean extends StsNumberFieldBean
{
	long minimum = -Long.MAX_VALUE;
	long maximum = Long.MAX_VALUE;
	long value = 0;
	int step = 1;
    byte oddState = ODD_STATE_NONE;

    static final byte ODD_STATE_NONE = 0;  // can be odd or even
    static final byte ODD_STATE_ALWAYS = 1;  // always odd
    static final byte ODD_STATE_NEVER = -1;  // always even

    boolean alwaysOdd = false;

    public StsLongFieldBean()
	{
	}

	public StsLongFieldBean(boolean editable)
	{
		this(editable, null);
	}

	public StsLongFieldBean(boolean editable, String label)
	{
		super(editable, label);
		layoutBean();
	}

	public StsLongFieldBean(Class c, String fieldName)
	{
        this.classInitialize(c, fieldName, true, null);
	}

	public StsLongFieldBean(Class c, String fieldName, String label)
	{
        this.classInitialize(c, fieldName, true, label);
	}

	public StsLongFieldBean(Class c, String fieldName, boolean editable, String label)
	{
        this.classInitialize(c, fieldName, editable, label);
	}

	public StsLongFieldBean(Class c, String fieldName, boolean editable)
	{
        this.classInitialize(c, fieldName, editable, null);
	}

	public StsLongFieldBean(Class c, String fieldName, long min, long max)
	{
        this.classInitialize(c, fieldName, min, max, true, null);
	}

	public StsLongFieldBean(Class c, String fieldName, long min, long max, String label)
	{
		this.classInitialize(c, fieldName, min, max, true, label);
	}

	public StsLongFieldBean(Class c, String fieldName, long min, long max, String label, boolean useArrows)
	{
		this.classInitialize(c, fieldName, min, max, true, label, useArrows);
	}

	public StsLongFieldBean(Class c, String fieldName, boolean editable, String label, boolean useArrows)
	{
		this.classInitialize(c, fieldName, editable, label, useArrows);
	}

	public StsLongFieldBean(Object beanObject, String fieldName, long min, long max, String label)
	{
		initialize(beanObject, fieldName, min, max, label);
	}

	public StsLongFieldBean(Object beanObject, String fieldName, boolean editable)
	{
		initialize(beanObject, fieldName, editable);
	}

	public StsLongFieldBean(Object beanObject, String fieldName, boolean editable, String label)
	{
		this.initialize(beanObject, fieldName, editable, label);
	}

	public StsLongFieldBean(Object beanObject, String fieldName, boolean editable, String label, boolean useArrows)
	{
		this.initialize(beanObject, fieldName, editable, label, useArrows);
	}

	public StsLongFieldBean(Object beanObject, String fieldName, long min, long max, String label, boolean useArrows)
	{
		super.initialize(beanObject, fieldName, true, label);
		setRange(min, max);
		layoutBean(useArrows);
    }

    public StsLongFieldBean(String labelText, long value)
    {
        setLabel(labelText);
        labelAndComponentLayout(label, textField);
        setValue(value);
        setEditable(false);
    }

	public void classInitialize(Class c, String fieldName, boolean editable, String label)
	{
		super.classInitialize(c, fieldName, editable, label);
        layoutBean();
    }

	public void classInitialize(Class c, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super.classInitialize(c, fieldName, editable, label);
		layoutBean(useArrows);
	}

    public void classInitialize(Class c, String fieldName, long min, long max, boolean editable, String label)
	{
		this.classInitialize(c, fieldName, min, max, editable, label, false);
	}

    public void classInitialize(Class c, String fieldName, long min, long max, boolean editable, String label, boolean useArrows)
	{
		super.classInitialize(c, fieldName, editable, label);
		setRange(min, max);
		layoutBean(useArrows);
//        addKeyListener(); // checks validity on each key stroke
	}

    public void initialize(Object beanObject, String fieldName, boolean editable)
	{
		this.initialize(beanObject, fieldName, editable, null);
        layoutBean();
    }

	public void initialize(Object beanObject, String fieldName, boolean editable, String label)
	{
		this.initialize(beanObject, fieldName, editable, label, false);
    }

	public void initialize(Object beanObject, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super.initialize(beanObject, fieldName, editable, label);
		layoutBean(useArrows);
    }

    public void initialize(Object beanObject, String fieldName, long min, long max, String label)
	{
		this.initialize(beanObject, fieldName, min, max, label, false);
	}

    public void initialize(Object beanObject, String fieldName, long min, long max, String label, boolean useArrows)
	{
		super.initialize(beanObject, fieldName, true, label);
		setRange(min, max);
		layoutBean(useArrows);
	}

	public StsLongFieldBean copy(Object beanObject)
	{
		StsLongFieldBean beanCopy = new StsLongFieldBean();
		if(editable)
			beanCopy.initialize(beanObject, fieldName, minimum, maximum, this.getLabelText(), useArrows);
		else
			beanCopy.initialize(beanObject, fieldName, false, getLabelText());
		return beanCopy;
	}

	private void setRange(long min, long max)
	{
		if(min <= max)
		{
			setMinimum(min);
			setMaximum(max);
		}
		else
		{
			setMinimum(max);
			setMaximum(min);
		}
	}

    public void setAlwaysOdd()
    {
        oddState = ODD_STATE_ALWAYS; 
        setStep(2);
    }
    public void setAlwaysEven()
    {
        oddState = ODD_STATE_NEVER;
        setStep(2);
    }

    // return true if value has changed
	public boolean setValue(String text)
	{
		try
		{
			long oldValue = value;
			value = Long.parseLong(text);
            return checkSetValue(value, text);
		}
		catch(Exception e)
		{
			StsToolkit.beep();
			textField.requestFocus();
			return false;
		}
	}

    private boolean checkSetValue(long value, String text)
    {
        boolean changed = false;
        long newValue = 0;
        if(oddState == ODD_STATE_NONE)
            newValue = value;
        if(oddState == ODD_STATE_ALWAYS)
            newValue = 2*(value/2) + 1;
        else if(oddState == ODD_STATE_NEVER)
            newValue = 2*(value/2);

        if(newValue != value)
        {
            value = newValue;
            text = Long.toString(value);
            changed = true;
        }
        setText(text);
        return changed;
    }

    public boolean setValue(long value)
	{
		long oldValue = this.value;
        if(value == oldValue) return false;
        this.value = value;
		super.setText(Long.toString(value));
		return value != oldValue;
	}

	public String toString()
	{
		return Long.toString(value);
	}

	public Object fromString(String string)
	{
		return Long.valueOf(string);
	}

	public void doSetValueObject(Object object)
	{
		long newValue = -1L;

		if(object instanceof String)
		{
			try { newValue = Long.parseLong((String)object); }
			catch(Exception e) { return; }
		}
		else if(object instanceof Long)
			newValue = ((Long)object).intValue();
		else if(object instanceof Integer)
			newValue = (long)((Integer)object).intValue();

//		if(value == newValue) return false;
		checkSetBoundedValue(newValue);
	}

	// return true if changed
	private boolean checkSetBoundedValue(long newValue)
	{
		if(newValue < minimum)
		{
			StsToolkit.beep();
			newValue = minimum;
		}
		else if(newValue > maximum)
		{
			StsToolkit.beep();
			newValue = maximum;
		}
		boolean changed = newValue != value;
		value = newValue;
        String text = Long.toString(value);
        checkSetValue(value, text);
		return changed;
	}

	public Object getValueObject()
	{
		return new Long(value);
	}
/*
	public Object getDefaultValueObject()
	{
		return new Long(0);
	}
*/
	public boolean increment()
	{
		long newValue = value + step;
		return checkSetBoundedValue(newValue);
	}

	public boolean decrement()
	{
		long newValue = value - step;
		return checkSetBoundedValue(newValue);
	}

	public long getMaximum() { return maximum; }
	public void setMaximum(long maximum) { this.maximum = maximum; }
	public long getMinimum() { return minimum; }
	public void setMinimum(long minimum) { this.minimum = minimum; }
	public int getStep() { return step; }
	public void setStep(int step) { this.step = step; }

	public long getValue()
	{
		return Long.parseLong(textField.getText());
	}

	public void setValueAndRange(long value, long min, long max)
	{
		if(max > min)
		{
			this.minimum = min;
			this.maximum = max;
		}
		else
		{
			this.maximum = min;
			this.minimum = max;
		}
		setValue(value);
	}


	public void setValueAndRangeFixStep(long value, long min, long max, long step)
	{
		setRange(min, max);
        setValue(value);
		setStep((int)step);
	}

    public boolean isValidChar(char c)
	{
		return Character.isDigit(c) || (c == KeyEvent.VK_BACK_SPACE) || (c == KeyEvent.VK_DELETE);
	}

	public boolean isValidString(String string)
	{
		if(string == null || string.length() <= 0) return true;

		try
		{
			Long.parseLong(string);
			return true;
		}
		catch(Exception e)
		{
			setText(Long.toString(value)); // restore value
			return false;
		}
	}
}
