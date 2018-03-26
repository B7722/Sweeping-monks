import QtQuick 2.7
import QtQuick.Controls 2.0
import QtQuick.VirtualKeyboard 2.1

ApplicationWindow {
    id: window
    visible: true
    width: 600
    height: 400
    title: qsTr("Stack")

    header: ToolBar {
        contentHeight: toolButton.implicitHeight

        ToolButton {
            id: toolButton
            text: stackView.depth > 1 ? "\u25C0" : "\u2630"
            font.pixelSize: Qt.application.font.pixelSize * 1.6
            onClicked: {
                if (stackView.depth > 1) {
                    stackView.pop()
                } else {
                    drawer.open()
                }
            }
        }

        Label {
            text: stackView.currentItem.title
            anchors.centerIn: parent
        }
    }

    Drawer {
        id: drawer
        width: window.width * 0.33
        height: window.height
        z: 100

        Column {
            anchors.fill: parent

            ItemDelegate {
                text: qsTr("Shape")
                width: parent.width
                onClicked: {
                    stackView.push("Page1Form.ui.qml")
                    drawer.close()
                }
            }
            ItemDelegate {
                text: qsTr("Enviroment")
                width: parent.width
                onClicked: {
                    stackView.push("Page2Form.ui.qml")
                    drawer.close()
                }
            }
            ItemDelegate {
                width: parent.width
                text: qsTr("Position")
                onClicked: {
                    stackView.push("Page3Form.ui.qml")
                    drawer.close()
                }
            }

            ItemDelegate {
                width: parent.width
                text: qsTr("Creator")
                onClicked: {
                    stackView.push("CreatorForm.ui.qml")
                    drawer.close()
                }
            }

        }
    }

    StackView {
        id: stackView
        initialItem: "HomeForm.ui.qml"
        anchors.fill: parent
    }

    InputPanel {
        id: inputPanel
        z: 99
        x: 0
        y: window.height
        width: window.width

        states: State {
            name: "visible"
            when: inputPanel.active
            PropertyChanges {
                target: inputPanel
                y: window.height - inputPanel.height
            }
        }
        transitions: Transition {
            from: ""
            to: "visible"
            reversible: true
            ParallelAnimation {
                NumberAnimation {
                    properties: "y"
                    duration: 250
                    easing.type: Easing.InOutQuad
                }
            }
        }
    }
}
