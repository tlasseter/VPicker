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

public class StsTextFieldBeanBeanInfo extends SimpleBeanInfo {
    private Class beanClass = StsTextFieldBean.class;
    private String iconColor16x16Filename;
    private String iconColor32x32Filename;
    private String iconMono16x16Filename;
    private String iconMono32x32Filename;

    public StsTextFieldBeanBeanInfo() {
    }
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            PropertyDescriptor _columns = new PropertyDescriptor("columns", beanClass, "getColumns", "setColumns");
            PropertyDescriptor _editable = new PropertyDescriptor("editable", beanClass, null, "setEditable");
            PropertyDescriptor _maximumSize = new PropertyDescriptor("maximumSize", beanClass, null, "setMaximumSize");
            PropertyDescriptor _minimumSize = new PropertyDescriptor("minimumSize", beanClass, null, "setMinimumSize");
            PropertyDescriptor _preferredSize = new PropertyDescriptor("preferredSize", beanClass, null, "setPreferredSize");
            PropertyDescriptor _text = new PropertyDescriptor("text", beanClass, "getText", "setText");
            PropertyDescriptor _valueFromParentObj = new PropertyDescriptor("valueFromParentObj", beanClass, null, "setValueFromParentObj");
            PropertyDescriptor[] pds = new PropertyDescriptor[] {
	            _columns,
	            _editable,
	            _maximumSize,
	            _minimumSize,
	            _preferredSize,
	            _text,
	            _valueFromParentObj};
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