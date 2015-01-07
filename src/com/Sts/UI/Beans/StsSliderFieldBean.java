package com.Sts.UI.Beans;

import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/** Constructs a JSlider which is connected to a given Class or instance of a Class.
 *  If a Class, the fieldName of its BoundedRangeModel must be provided;
 *  if an instance of a class, the BoundedRangeModel itself is passed in as an argument.
 *  Both constructors require a changedMethodName which is the name of a method in the Class
 *  or instance of a Class which is called without arguments when the rangeModel state is
 *  changed.  The called method should then examine the rangeModel and take appropriate actions.
 */
public class StsSliderFieldBean extends StsFieldBean implements ChangeListener
{
	private JSlider slider = new JSlider();
	private StsIntFieldBean valueBean;
    StsMethod changedMethod;

    public StsSliderFieldBean()
    {
    }

    /**
     * Constructor for a Class which has a BoundedRangeModel as a field and a rangeModel changed method.
     * Subsequently, the beanObject which is an instance of this Class must be defined.
     * @param c Class of instance containing a BoundedRangeModel and changedModel methodName.
     * @param rangeModelFieldname fieldname of the BoundedRangeModel for this class; when beanObject is set, rangeModel is set in JSlider.
     * @param changedMethodName called whenever the fireStateChanged method is called by JSlider.
     */
    public StsSliderFieldBean(Class c, String rangeModelFieldname, String changedMethodName)
    {
		super();
        classInitialize(c, rangeModelFieldname, false, null);
        changedMethod = new StsMethod(c, changedMethodName);
        initialize();
    }

    /**
     * Constructor for an instance of a Class which has a BoundedRangeModel as a field and a rangeModel changed method.
     * @param instance instance containing a BoundedRangeModel and changedModel methodName.
     * @param rangeModel BoundedRangeModel for an instance
     * @param changedMethodName called whenever the fireStateChanged method is called by JSlider.
     */
    public StsSliderFieldBean(Object instance, BoundedRangeModel rangeModel, String changedMethodName)
    {
		super();
        beanObject = instance;
        changedMethod = new StsMethod(instance, changedMethodName);
        initialize();
        setRangeModel(rangeModel);
    }

    public void initialize()
    {
    	valueBean = new StsIntFieldBean();
    	valueBean.initialize(this, "Value", false);
        slider.setPaintTicks(true);
        slider.setPaintLabels(true);
        slider.addChangeListener(this);
        gbc.fill = gbc.HORIZONTAL;
        gbc.weightx = 1.0;
        addToRow(slider);
    }

    public void setRangeModel(BoundedRangeModel rangeModel)
    {
        slider.setModel(rangeModel);
        int min = rangeModel.getMinimum();
        int max = rangeModel.getMaximum();
        double[] scale = StsMath.niceScaleBold(min, max, 10, true);
        int minorTickSpacing = (int)scale[2];
        int majorTickSpacing = (int)scale[3];
        slider.setMajorTickSpacing(majorTickSpacing);
        slider.setMinorTickSpacing(minorTickSpacing);
    }

    public void setVertical() { slider.setOrientation(JSlider.VERTICAL); }
    public void setHorizontal() { slider.setOrientation(JSlider.HORIZONTAL); }
    public void setOrientation(int orientation) { slider.setOrientation(orientation); }
    public void setMajorTickSpacing(int spacing) { slider.setMajorTickSpacing(spacing); }
    public void setMinorTickSpacing(int spacing) { slider.setMinorTickSpacing(spacing); }
    public void setPaintTicks(boolean paint) { slider.setPaintTicks(paint); }
    public void setPaintLabels(boolean paint) { slider.setPaintLabels(paint); }
    public void setInverted(boolean inverted) { slider.setInverted(inverted); }
    public void setShowValue(boolean showValue)
    {
    	//gbc.fill = gbc.HORIZONTAL;   	
    	addEndRow(valueBean);
    	valueBean.setValue(String.valueOf(slider.getValue()));
    }

    public Component[] getBeanComponents() { return new Component[] { slider }; }

	public String toString() { return StsParameters.NONE_STRING; }
	public Object fromString(String string) { return null; }

    public Object getValueObject() { return getModel(); }
    public BoundedRangeModel getModel() { return slider.getModel(); }

	public float getValue() { return (float)getModel().getValue(); }

    public void doSetValueObject(Object object)
    {
        if(!(object instanceof BoundedRangeModel))
        {
            super.outputSetValueObjectException(object, "BoundedRangeModel");
            return;
        }
        BoundedRangeModel rangeModel = (BoundedRangeModel)object;
        setRangeModel(rangeModel);
        slider.setModel(rangeModel);
    }

    public void stateChanged(ChangeEvent e)
    {
        if(beanObject == null)
            changedMethod.invokeStaticMethod(null);
        else
            changedMethod.invokeInstanceMethod(beanObject, new Object[0]);
    	valueBean.setValue(String.valueOf(slider.getValue()));
    }

    public void actionPerformed(ActionEvent e)
    {
        System.out.println("slider actionPerformed called.");
    }

    static public void main(String[] args)
    {
        SliderTest sliderTest = new SliderTest(0, 10000, 2500);
        StsSliderFieldBean sliderBean = new StsSliderFieldBean(sliderTest, sliderTest.rangeModel,  "rangeModelChanged");
        sliderBean.setOrientation(JSlider.HORIZONTAL);
//        StsSliderFieldBean sliderBean = new StsSliderFieldBean(SliderTest.class, "rangeModel",  "rangeModelChanged");
//        sliderBean.setBeanObject(sliderTest);
        StsToolkit.createDialog(sliderBean);
    }
}
