#include <avr/io.h>
#include <avr/pgmspace.h>
#include <util/delay.h>
#include "usb_rawhid.h"

#define LED_ON      (PORTE |= _BV(PE6))
#define LED_OFF     (PORTE &= ~_BV(PE6))
#define LED_TOGGLE  (PORTE ^= _BV(PE6))

uint8_t rec_buf[8];
uint8_t send_buf[] = "USB on an Atmega!!";

int main()
{
    DDRE = 0xff;
    usb_init();
    while(!usb_configured());
    
    uint8_t r = 1;
    while(1)
    {
        _delay_ms(5000);
        //r = usb_rawhid_recv(rec_buf, 0);
        if( r > 0)
        {
            LED_TOGGLE;
            usb_rawhid_send(send_buf, 10);
        }
    }
    return 0;
}
