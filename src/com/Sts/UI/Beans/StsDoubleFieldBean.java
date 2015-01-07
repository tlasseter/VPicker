package com.Sts.UI.Beans;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

import com.Sts.Utilities.*;

import java.text.*;

public class StsDoubleFieldBean extends StsNumberFieldBean
{
    double minimum = -Double.MAX_VALUE;
    double maximum = Double.MAX_VALUE;
    double value = minimum;
    double step = 1.0;
	boolean stepFixed = false;
    DecimalFormat labelFormat = new DecimalFormat();

    public StsDoubleFieldBean()
    {
    }

    public StsDoubleFieldBean(boolean editable)
    {
        this(editable, null);
    }

    public StsDoubleFieldBean(boolean editable, String label)
    {
        super(editable, label);
		layoutBean();
    }

    public StsDoubleFieldBean(Class c, String fieldName)
    {
        this.classInitialize(c, fieldName, true, null);
    }

    public StsDoubleFieldBean(Class c, String fieldName, String label)
    {
        this.classInitialize(c, fieldName, true, label);
    }

    public StsDoubleFieldBean(Class c, String fieldName, boolean editable, String label)
    {
        this.classInitialize(c, fieldName, editable, label);
    }

    public StsDoubleFieldBean(Class c, String fieldName, boolean editable)
    {
        this.classInitialize(c, fieldName, editable, null);
    }

    public StsDoubleFieldBean(Class c, String fieldName, double min, double max)
    {
        this.classInitialize(c, fieldName, min, max, true, null);
    }

    public StsDoubleFieldBean(Class c, String fieldName, double value)
    {
        this.classInitialize(c, fieldName, true, null);
    }

    public StsDoubleFieldBean(Class c, String fieldName, double min, double max, String label)
    {
        this.classInitialize(c, fieldName, min, max, true, label);
    }

	public StsDoubleFieldBean(Class c, String fieldName, boolean editable, String label, boolean useArrows)
	{
		this.classInitialize(c, fieldName, editable, label, useArrows);
	}

    public StsDoubleFieldBean(Object beanObject, String fieldName, String label)
    {
        this.initialize(beanObject, fieldName, true, label);
    }

    public StsDoubleFieldBean(Object beanObject, String fieldName)
    {
        this(beanObject, fieldName, null);
    }

	public StsDoubleFieldBean(Object beanObject, String fieldName, boolean editable, String label)
	{
		this(beanObject, fieldName, editable, label, false);
	}

	public StsDoubleFieldBean(Object beanObject, String fieldName, boolean editable, String label, boolean useArrows)
	{
		this.initialize(beanObject, fieldName, editable, label, useArrows);
   }

	public StsDoubleFieldBean(Object beanObject, String fieldName, double min, double max, String label)
	{
		this.initialize(beanObject, fieldName, min, max, label);
    }

	public StsDoubleFieldBean(Object beanObject, String fieldName, double min, double max, String label, boolean useArrows)
	{
		this.initialize(beanObject, fieldName, min, max, label, useArrows);
   }

	public void classInitialize(Class c, String fieldName, boolean editable, String label)
	{
		this.classInitialize(c, fieldName, editable, label, false);
    }

	public void classInitialize(Class c, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super.classInitialize(c, fieldName, editable, label);
        layoutBean(useArrows);
    }

	public void classInitialize(Class c, String fieldName, double min, double max, String label)
    {
        this.classInitialize(c, fieldName, min, max, true, label);
    }

    public void classInitialize(Class c, String fieldName, double min, double max, boolean editable, String label)
    {
        this.classInitialize(c, fieldName, min, max, editable, label, false);
    }

	public void classInitialize(Class c, String fieldName, double min, double max, boolean editable, String label, boolean useArrows)
	{
		super.classInitialize(c, fieldName, editable, label);
		setRange(min, max);
		layoutBean(useArrows);
    }

	public void initialize(Object beanObject, String fieldName, boolean editable)
	{
		this.initialize(beanObject, fieldName, editable, null);
        layoutBean();
    }

	public void initialize(Object beanObject, String fieldName, boolean editable, String label)
	{
		super.initialize(beanObject, fieldName, editable, label);
		layoutBean();
    }

	public void initialize(Object beanObject, String fieldName, boolean editable, String label, boolean useArrows)
	{
		super.initialize(beanObject, fieldName, editable, label);
		layoutBean(useArrows);
    }

    public void initialize(Object beanObject, String fieldName, double min, double max, String label)
    {
        this.initialize(beanObject, fieldName, min, max, label, false);
    }

	public void initialize(Object beanObject, String fieldName, double min, double max, String label, boolean useArrows)
	{
		super.initialize(beanObject, fieldName, true, label);
		setRange(min, max);
		layoutBean(useArrows);
    }

	public StsDoubleFieldBean copy(Object beanObject)
	{
		StsDoubleFieldBean beanCopy = new StsDoubleFieldBean();
		if(editable)
			beanCopy.initialize(beanObject, fieldName, minimum, maximum, this.getLabelText(), useArrows);
		else
			beanCopy.initialize(beanObject, fieldName, false, getLabelText());
		return beanCopy;
	}

    public void setRange(double min, double max)
    {
		if(min <= max)
		{
			setMinimum(min);
			setMaximum(max);
			step = StsMath.niceNumber((max - min)/100,false);
			fixStep(step, false);
		}
		else
		{
			setMinimum(max);
			setMaximum(min);
			step = StsMath.niceNumber((min - max)/100, true);
			fixStep(step, false);
		}
	}

	public void setValueAndRange(double value, double min, double max)
	{
		setRange(min, max);
		setValue(value);
	}

	public void setValueAndRangeFixStep(double value, double min, double max, double step)
	{
		setRange(min, max);
        setValue(value);
		fixStep(step, true);
	}

	public void setRangeFixStep(double min, double max, double step)
	{
		setRange(min, max);
		fixStep(step, true);
	}

    public void fixStep(double step)
	{
		fixStep(step, true);
	}

	protected void fixStep(double step, boolean stepFixed)
	{
		this.stepFixed = stepFixed;
		this.step = Math.abs(step);
	/*
		Math.log(step);
		double pow = StsMath.log10(step);
		int powi = (int)Math.round(pow);
		if(pow < 0)
            labelFormat.setMaximumFractionDigits(-powi);
		step = Math.pow(10, powi);
	*/
	}

    public boolean setValue(double value)
    {
		if(this.value == value) return false;
        return checkSetBoundedValue(value);
    }

    public boolean setValueObject(Object object)
    {
        double newValue;

        if(object instanceof String)
        {
            try { newValue = labelFormat.parse((String)object).doubleValue(); }
            catch(Exception e) { return false; }
        }
        else if(object instanceof Double)
            newValue = ((Double)object).doubleValue();
		else if(object instanceof Float)
			newValue = ((Float)object).doubleValue();
        else
            return false;

//        if(value == newValue) return false;
        return checkSetBoundedValue(newValue);
    }

    public void setSignificantDigits(int digits)
    {
		labelFormat.setMinimumIntegerDigits(1);
		labelFormat.setMaximumIntegerDigits(4);
		labelFormat.setMinimumFractionDigits(4);
    }

    // return true if changed
    private boolean checkSetBoundedValue(double newValue)
    {
        if((float)newValue < (float)minimum)
        {
            StsToolkit.beep();
            newValue = minimum;
        }
        else if((float)newValue > (float)maximum)
        {
            StsToolkit.beep();
            newValue = maximum;
        }

	    if (stepFixed)
		{
			newValue = this.computeIncrementalValue(newValue);
//			setText(Double.toString(newValue));
			setText(formatNumber(newValue));
		}
		else
//			setText(Double.toString(newValue));
            setText(formatNumber(newValue));

		value = newValue;
		return true;
    }

    public void setFormat(DecimalFormat format)
    {
        labelFormat = format;
    }

    public String formatNumber(double number)
	{
		return formatNumber(number, 10);
	}

	public String formatNumber(double number, int nSigDigits)
	{
		if(number == 0) return new String("0.0");

		int nIntegerDigits = (int)StsMath.log10(Math.abs(number)) + 1;
		boolean exponentialPattern = nIntegerDigits > nSigDigits || nIntegerDigits < -nSigDigits;
		DecimalFormat format = new DecimalFormat();

		if(exponentialPattern)
		{
			format.applyPattern("0.E0");
			format.setMinimumIntegerDigits(1);
			format.setMaximumFractionDigits(nSigDigits - 1);
		}
		else if(nIntegerDigits > 0)
		{
			format.setMinimumFractionDigits(0);
			format.setMaximumFractionDigits(StsMath.minMax(nSigDigits - nIntegerDigits, 0, 2));
		}
		else
		{
			format.setMinimumFractionDigits(0);
			format.setMaximumFractionDigits(nSigDigits - nIntegerDigits);
		}
		String numberString = format.format(number);
		int nChars = numberString.length();
		if(nChars > nColumns && !exponentialPattern)
		{
//			System.out.println("string " + numberString + " larger than " + nColumns + " columns. Revising...");
			format.applyPattern("0.E0");
			format.setMinimumIntegerDigits(1);
			format.setMaximumFractionDigits(Math.max(0, nColumns - 3));
		}
		//System.out.println(format.toPattern() + " number: " + number + " string: " + numberString);
		return numberString;
	}

    public Object getValueObject()
    {
        return new Double(value);
    }

	public String toString()
	{
		return Double.toString(value);
	}

	public Object fromString(String string)
	{
        return Double.valueOf(string);
	}
/*
	public Object getDefaultValueObject()
	{
		return new Double(0);
	}
*/
    public boolean increment()
    {
        double newValue = value + step;
        return checkSetBoundedValue(newValue);
    }

	private double computeIncrementalValue(double value)
	{
		double steps;
		if (minimum != -Double.MAX_VALUE)
			steps = (value - minimum) / step;
		else
			steps = value / step;

		long isteps = Math.round(steps);
//		double remainder = Math.abs(steps - isteps);
//            if(remainder > 1.e-5) StsToolkit.beep();
		if (minimum != -Double.MAX_VALUE)
			value = isteps * step + minimum;
		else
			value = isteps * step;
		return value;
	}

    public boolean decrement()
    {
        double newValue = value - step;
        return checkSetBoundedValue(newValue);
    }

    public double getMaximum() { return maximum; }
    public void setMaximum(double maximum) { this.maximum = maximum; }
    public double getMinimum() { return minimum; }
    public void setMinimum(double minimum) { this.minimum = minimum; }

    public double getValue()
    {
        return value;
    }

	/** String has been directly typed into textField */
    public boolean isValidString(String string)
    {
        if(string == null || string.length() <= 0) return true;
		try
		{
			double newValue = labelFormat.parse(string).doubleValue();
//			double newValue = labelFormat.pDouble.parseDouble(string);
//			setText(string);
//			value = newValue;
			return true;
		}
		catch(Exception e)
		{
			setValue(value);
			return false;
		}
    }
}
