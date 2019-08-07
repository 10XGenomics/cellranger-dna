<?php
// DO NOT EDIT! Generated by Protobuf-PHP protoc plugin 1.0
// Source: tabletmanagerdata.proto

namespace Vitess\Proto\Tabletmanagerdata {

  class StopSlaveMinimumRequest extends \DrSlump\Protobuf\Message {

    /**  @var string */
    public $position = null;
    
    /**  @var int */
    public $wait_timeout = null;
    

    /** @var \Closure[] */
    protected static $__extensions = array();

    public static function descriptor()
    {
      $descriptor = new \DrSlump\Protobuf\Descriptor(__CLASS__, 'tabletmanagerdata.StopSlaveMinimumRequest');

      // OPTIONAL STRING position = 1
      $f = new \DrSlump\Protobuf\Field();
      $f->number    = 1;
      $f->name      = "position";
      $f->type      = \DrSlump\Protobuf::TYPE_STRING;
      $f->rule      = \DrSlump\Protobuf::RULE_OPTIONAL;
      $descriptor->addField($f);

      // OPTIONAL INT64 wait_timeout = 2
      $f = new \DrSlump\Protobuf\Field();
      $f->number    = 2;
      $f->name      = "wait_timeout";
      $f->type      = \DrSlump\Protobuf::TYPE_INT64;
      $f->rule      = \DrSlump\Protobuf::RULE_OPTIONAL;
      $descriptor->addField($f);

      foreach (self::$__extensions as $cb) {
        $descriptor->addField($cb(), true);
      }

      return $descriptor;
    }

    /**
     * Check if <position> has a value
     *
     * @return boolean
     */
    public function hasPosition(){
      return $this->_has(1);
    }
    
    /**
     * Clear <position> value
     *
     * @return \Vitess\Proto\Tabletmanagerdata\StopSlaveMinimumRequest
     */
    public function clearPosition(){
      return $this->_clear(1);
    }
    
    /**
     * Get <position> value
     *
     * @return string
     */
    public function getPosition(){
      return $this->_get(1);
    }
    
    /**
     * Set <position> value
     *
     * @param string $value
     * @return \Vitess\Proto\Tabletmanagerdata\StopSlaveMinimumRequest
     */
    public function setPosition( $value){
      return $this->_set(1, $value);
    }
    
    /**
     * Check if <wait_timeout> has a value
     *
     * @return boolean
     */
    public function hasWaitTimeout(){
      return $this->_has(2);
    }
    
    /**
     * Clear <wait_timeout> value
     *
     * @return \Vitess\Proto\Tabletmanagerdata\StopSlaveMinimumRequest
     */
    public function clearWaitTimeout(){
      return $this->_clear(2);
    }
    
    /**
     * Get <wait_timeout> value
     *
     * @return int
     */
    public function getWaitTimeout(){
      return $this->_get(2);
    }
    
    /**
     * Set <wait_timeout> value
     *
     * @param int $value
     * @return \Vitess\Proto\Tabletmanagerdata\StopSlaveMinimumRequest
     */
    public function setWaitTimeout( $value){
      return $this->_set(2, $value);
    }
  }
}

