<?php
// DO NOT EDIT! Generated by Protobuf-PHP protoc plugin 1.0
// Source: tabletmanagerdata.proto

namespace Vitess\Proto\Tabletmanagerdata {

  class SchemaDefinition extends \DrSlump\Protobuf\Message {

    /**  @var string */
    public $database_schema = null;
    
    /**  @var \Vitess\Proto\Tabletmanagerdata\TableDefinition[]  */
    public $table_definitions = array();
    
    /**  @var string */
    public $version = null;
    

    /** @var \Closure[] */
    protected static $__extensions = array();

    public static function descriptor()
    {
      $descriptor = new \DrSlump\Protobuf\Descriptor(__CLASS__, 'tabletmanagerdata.SchemaDefinition');

      // OPTIONAL STRING database_schema = 1
      $f = new \DrSlump\Protobuf\Field();
      $f->number    = 1;
      $f->name      = "database_schema";
      $f->type      = \DrSlump\Protobuf::TYPE_STRING;
      $f->rule      = \DrSlump\Protobuf::RULE_OPTIONAL;
      $descriptor->addField($f);

      // REPEATED MESSAGE table_definitions = 2
      $f = new \DrSlump\Protobuf\Field();
      $f->number    = 2;
      $f->name      = "table_definitions";
      $f->type      = \DrSlump\Protobuf::TYPE_MESSAGE;
      $f->rule      = \DrSlump\Protobuf::RULE_REPEATED;
      $f->reference = '\Vitess\Proto\Tabletmanagerdata\TableDefinition';
      $descriptor->addField($f);

      // OPTIONAL STRING version = 3
      $f = new \DrSlump\Protobuf\Field();
      $f->number    = 3;
      $f->name      = "version";
      $f->type      = \DrSlump\Protobuf::TYPE_STRING;
      $f->rule      = \DrSlump\Protobuf::RULE_OPTIONAL;
      $descriptor->addField($f);

      foreach (self::$__extensions as $cb) {
        $descriptor->addField($cb(), true);
      }

      return $descriptor;
    }

    /**
     * Check if <database_schema> has a value
     *
     * @return boolean
     */
    public function hasDatabaseSchema(){
      return $this->_has(1);
    }
    
    /**
     * Clear <database_schema> value
     *
     * @return \Vitess\Proto\Tabletmanagerdata\SchemaDefinition
     */
    public function clearDatabaseSchema(){
      return $this->_clear(1);
    }
    
    /**
     * Get <database_schema> value
     *
     * @return string
     */
    public function getDatabaseSchema(){
      return $this->_get(1);
    }
    
    /**
     * Set <database_schema> value
     *
     * @param string $value
     * @return \Vitess\Proto\Tabletmanagerdata\SchemaDefinition
     */
    public function setDatabaseSchema( $value){
      return $this->_set(1, $value);
    }
    
    /**
     * Check if <table_definitions> has a value
     *
     * @return boolean
     */
    public function hasTableDefinitions(){
      return $this->_has(2);
    }
    
    /**
     * Clear <table_definitions> value
     *
     * @return \Vitess\Proto\Tabletmanagerdata\SchemaDefinition
     */
    public function clearTableDefinitions(){
      return $this->_clear(2);
    }
    
    /**
     * Get <table_definitions> value
     *
     * @param int $idx
     * @return \Vitess\Proto\Tabletmanagerdata\TableDefinition
     */
    public function getTableDefinitions($idx = NULL){
      return $this->_get(2, $idx);
    }
    
    /**
     * Set <table_definitions> value
     *
     * @param \Vitess\Proto\Tabletmanagerdata\TableDefinition $value
     * @return \Vitess\Proto\Tabletmanagerdata\SchemaDefinition
     */
    public function setTableDefinitions(\Vitess\Proto\Tabletmanagerdata\TableDefinition $value, $idx = NULL){
      return $this->_set(2, $value, $idx);
    }
    
    /**
     * Get all elements of <table_definitions>
     *
     * @return \Vitess\Proto\Tabletmanagerdata\TableDefinition[]
     */
    public function getTableDefinitionsList(){
     return $this->_get(2);
    }
    
    /**
     * Add a new element to <table_definitions>
     *
     * @param \Vitess\Proto\Tabletmanagerdata\TableDefinition $value
     * @return \Vitess\Proto\Tabletmanagerdata\SchemaDefinition
     */
    public function addTableDefinitions(\Vitess\Proto\Tabletmanagerdata\TableDefinition $value){
     return $this->_add(2, $value);
    }
    
    /**
     * Check if <version> has a value
     *
     * @return boolean
     */
    public function hasVersion(){
      return $this->_has(3);
    }
    
    /**
     * Clear <version> value
     *
     * @return \Vitess\Proto\Tabletmanagerdata\SchemaDefinition
     */
    public function clearVersion(){
      return $this->_clear(3);
    }
    
    /**
     * Get <version> value
     *
     * @return string
     */
    public function getVersion(){
      return $this->_get(3);
    }
    
    /**
     * Set <version> value
     *
     * @param string $value
     * @return \Vitess\Proto\Tabletmanagerdata\SchemaDefinition
     */
    public function setVersion( $value){
      return $this->_set(3, $value);
    }
  }
}

