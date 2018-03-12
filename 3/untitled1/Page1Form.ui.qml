import QtQuick 2.7
import QtQuick.Controls 2.0

Page {
    width: 600
    height:400
    background: black
    title: qsTr("Shape")

    Button {
        id: button
        x: 452
        y: 299
        width: 102
        height: 32
        text: qsTr("Cancel")
    }

    Label {
        id: label2
        x: 56
        y: 60
        width: 122
        height: 32
        text: qsTr("Parameters")
    }

    Label {
        id: label3
        x: 56
        y: 113
        width: 122
        height: 32
        text: qsTr("Degree")
    }

    ComboBox {
        id: comboBox
        x: 212
        y: 8
    }

    ComboBox {
        //model: [ "Banana", "Apple", "Coconut" ]
        id: comboBox1
        x: 212
        y: 60
    }

    TextField {
        id: textField
        x: 435
        y: 60
        width:119
        height: 32
    }

    Label {
        id: label4
        x: 381
        y: 66
        width: 64
        height: 32
        text: qsTr("Value")
    }

    Button {
        id: button1
        x: 320
        y: 299
        width: 102
        height: 32
        text: qsTr("Set")
    }

    Label {
        id: label5
        x: 56
        y: 0
        width: 122
        height: 32
        text: "Shape"
        verticalAlignment: Text.AlignVCenter
        horizontalAlignment: Text.AlignLeft
        fontSizeMode: Text.FixedSize
    }

    Label {
        id: label6
        x: 56
        y: 166
        width: 122
        height: 32
        text: "Camera"
        verticalAlignment: Text.AlignVCenter
        horizontalAlignment: Text.AlignLeft
        fontSizeMode: Text.FixedSize
    }

    ComboBox {
        id: comboBox2
        x: 212
        y: 166
    }

    TextField {
        id: textField1
        x: 213
        y: 113
        width: 119
        height: 32
    }
}
