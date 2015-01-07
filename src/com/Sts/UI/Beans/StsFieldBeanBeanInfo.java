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

public class StsFieldBeanBeanInfo extends SimpleBeanInfo {
    private Class beanClass = StsFieldBean.class;
    private String iconColor16x16Filename;
    private String iconColor32x32Filename;
    private String iconMono16x16Filename;
    private String iconMono32x32Filename;

    public StsFieldBeanBeanInfo() {
    }
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            PropertyDescriptor _fieldName = new PropertyDescriptor("fieldName", beanClass, "getFieldName", "setFieldName");
            PropertyDescriptor _getMethod = new PropertyDescriptor("getMethod", beanClass, "getGetMethod", "setGetMethod");
            PropertyDescriptor _setMethod = new PropertyDescriptor("setMethod", beanClass, "getSetMethod", "setSetMethod");
            PropertyDescriptor _text = new PropertyDescriptor("text", beanClass, "getText", null);
            PropertyDescriptor _updateMethod = new PropertyDescriptor("updateMethod", beanClass, "getUpdateMethod", null);
            PropertyDescriptor _valueFromParentObj = new PropertyDescriptor("valueFromParentObj", beanClass, null, "setValueFromParentObj");
            PropertyDescriptor _valueObject = new PropertyDescriptor("valueObject", beanClass, "getValueObject", "setValueObject");
            PropertyDescriptor[] pds = new PropertyDescriptor[] {
	            _fieldName,
	            _getMethod,
	            _setMethod,
	            _text,
	            _updateMethod,
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