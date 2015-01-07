package com.Sts.UI.Beans;

import java.beans.*;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

public class StsFloatFieldBeanBeanInfo extends SimpleBeanInfo {
    private Class beanClass = StsFloatFieldBean.class;
    private String iconColor16x16Filename;
    private String iconColor32x32Filename = "StsFloatFieldBean.gif";
    private String iconMono16x16Filename;
    private String iconMono32x32Filename;

    public StsFloatFieldBeanBeanInfo() {
    }
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            PropertyDescriptor _maximum = new PropertyDescriptor("maximum", beanClass, "getMaximum", "setMaximum");
            PropertyDescriptor _minimum = new PropertyDescriptor("minimum", beanClass, "getMinimum", "setMinimum");
            PropertyDescriptor _step = new PropertyDescriptor("step", beanClass, "getStep", "setStep");
            PropertyDescriptor _value = new PropertyDescriptor("value", beanClass, "getValue", "setValue");
            PropertyDescriptor _valueFromParentObj = new PropertyDescriptor("valueFromParentObj", beanClass, null, "setValueFromParentObj");
            PropertyDescriptor _valueObject = new PropertyDescriptor("valueObject", beanClass, "getValueObject", "setValueObject");
            PropertyDescriptor[] pds = new PropertyDescriptor[] {
	            _maximum,
	            _minimum,
	            _step,
	            _value,
	            _valueFromParentObj,
	            _valueObject};
            return pds;
        }
        catch(IntrospectionException ex) {
            ex.printStackTrace();
            return null;
        }
    }
    public java.awt.Image getIcon(int iconKind) {
        switch (iconKind) {
            case BeanInfo.ICON_COLOR_16x16:
              return iconColor16x16Filename != null ? loadImage(iconColor16x16Filename) : null;
            case BeanInfo.ICON_COLOR_32x32:
              return iconColor32x32Filename != null ? loadImage(iconColor32x32Filename) : null;
            case BeanInfo.ICON_MONO_16x16:
              return iconMono16x16Filename != null ? loadImage(iconMono16x16Filename) : null;
            case BeanInfo.ICON_MONO_32x32:
              return iconMono32x32Filename != null ? loadImage(iconMono32x32Filename) : null;
        }
        return null;
    }
    public BeanInfo[] getAdditionalBeanInfo() {
        Class superclass = beanClass.getSuperclass();
        try {
            BeanInfo superBeanInfo = Introspector.getBeanInfo(superclass);
            return new BeanInfo[] { superBeanInfo };
        }
        catch(IntrospectionException ex) {
            ex.printStackTrace();
            return null;
        }
    }
}
